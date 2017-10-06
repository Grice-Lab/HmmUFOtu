/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/*
 * BandedHMMP7.h
 *
 *  Created on: May 13, 2015
 *      Author: zhengqi
 */

#ifndef BANDEDHMMP7_H_
#define BANDEDHMMP7_H_
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <Eigen/Dense>
#include <limits>
#include <cmath>
#include <cstdio>
#include <climits>
#include <stdint.h> /* for fixed size integers */
#include <iostream>
#include "HmmUFOtuConst.h"
#include "HmmUFOtuDef.h"
#include "AlphabetFactory.h"
#include "BandedHMMP7Bg.h"
#include "BandedHMMP7Prior.h"
#include "StringUtils.h"
#include "PrimarySeq.h"
#include "DigitalSeq.h"
#include "MSA.h"
#include "CSLoc.h"
#include "RootFinder.h"

namespace EGriceLab {
using std::string;
using std::istream;
using std::ostream;
using std::vector;
using std::map;
using std::deque;
using Eigen::Matrix3d;
using Eigen::Matrix4Xd;
using Eigen::Matrix4d;
using Math::RootFinder;

/**
 * Banded plan7 HMM for 16S rRNA profile alignment
 * Similar to HMMER, the global profile includes N, B, M, I, D, E, C states
 * but not the J (joining) state, so repeated and multiple mapping of same region not allowed
 */
class BandedHMMP7 {
public:
	/* constructors */
	/**
	 * Default constructor, do zero initiation
	 */
	BandedHMMP7();
	/**
	 * Construct a BandedHMMP7 with given length and alphabet
	 */
	BandedHMMP7(const string& name, int K, const DegenAlphabet* abc);

	/**
	 * Construct a BandedHMMP7 with given version, length and alphabet
	 */
	BandedHMMP7(const string& name, const string& version, int K, const DegenAlphabet* abc);

	/* nested enums and types */
	/* enum members of all P7 states
	 * M: match
	 * I: insertion
	 * D: deletion
	 * N: N or 5'
	 * C: C or 3'
	 * B: begin state
	 * E: end state
	 */
	enum p7_state { M, I, D, N, C, B, E, P /* non-exist phantom state */ };

	/** align mode relative to the read */
	enum align_mode {
		GLOBAL,
		LOCAL,
		NGCL /* N' global C' local */,
		CGNL /* C' global N' local */
	};

	/** padding mode for filling non profile CS positions */
	enum padding_mode {
		LEFT,
		RIGHT,
		MIDDLE,
		JUSTIFIED
	};

	/* forward declaration of nested classes and alias */
	struct ViterbiScores; /* struct storing the ViterbiScores used during the Viterbi algorithm */
	struct ViterbiAlignPath; /* a known viterbi align path */
	struct ViterbiAlignTrace; /* the final backtrace of the Viterbi algorithm */

	typedef ViterbiScores VScore;
	typedef ViterbiAlignPath VPath;
	typedef ViterbiAlignTrace VTrace;

	/**
	 * A nested class storing public accessible Viterbi Scoring matrices
	 */
	struct ViterbiScores {
		/* constructors*/
		/** default constructor, do nothing */
		ViterbiScores() : K(0), L(0) {  }

		/** construct a VScore with given sizes */
		explicit ViterbiScores(int K, int L) : K(K), L(L)
		{
			reset();
		}

		/* member fields */
		const int K; /* fixed profile size */
		int L; /* current seq length */
		/* Viterbi cost matrices */
		MatrixXd DP_M;  /* (L+1) * (K+1) cost matrix of the best path matching the subsequence X1..i
							to the profile submodel up to the column j, ending with xi being emitted by Mj*/
		MatrixXd DP_I;  /* (L+1) * (K+1) cost matrix of the best path matching the subsequence Xi..i
							to the profile submodel up to the ending in xi being emitted by Ij */
		MatrixXd DP_D;  /* (L+1) * (K+1) cost matrix of the best path ending in Dj, and xi being the last character emitted before Dj).*/

		MatrixXd S;     /* (L+1) * (K+2) score matrix storing the whole ViterbiScores with the last column the cost exiting from Dk status */

		/* member methods */
		void reset(int L) {
			this->L = L;
			reset();
		}

		/** reset all matrix values according the current seq length */
		void reset() {
			DP_M.resize(L + 1, K + 1);
			DP_I.resize(L + 1, K + 1);
			DP_D.resize(L + 1, K + 1);
			S.resize(L + 1, K + 2);

			DP_M.setConstant(inf);
			DP_I.setConstant(inf);
			DP_D.setConstant(inf);
			S.setConstant(inf);
		}
	};

	struct ViterbiAlignPath {
		/* constructors */
		/** default constructor, do nothing */
		ViterbiAlignPath() { }

		/** construct a VPath from given information */
		ViterbiAlignPath(int start, int end, int from, int to, int nIns, int nDel) :
			start(start), end(end), from(from), to(to), nIns(nIns), nDel(nDel)
		{  }

		/** member methods */
		bool isValid() const {
			return start > 0 && start <= end && from > 0 && from <= to && nIns >= 0 && nDel >= 0;
		}

		/* member fields */
		int start, end; /* 1-based position on profile */
		int from, to;   /* 1-based position on seq */
		int nIns;       /* known number of insertions */
		int nDel;       /* known number of deletions */
	};

	struct ViterbiAlignTrace {
		/* constructors */
		/** default constructor */
		ViterbiAlignTrace() :
			minScore(inf), alnStart(0), alnEnd(0), alnFrom(0), alnTo(0)
		{  }

		/** construct a VTrace using a given initial information */
		ViterbiAlignTrace(double minScore, int alnStart, int alnEnd, int alnFrom, int alnTo, string alnTrace) :
			minScore(minScore), alnStart(alnStart), alnEnd(alnEnd), alnFrom(alnFrom), alnTo(alnTo), alnTrace(alnTrace)
		{  }

		/* member fields */
		double minScore;
		int alnStart, alnEnd; // final align start and end (1-based)
		int alnFrom, alnTo; // final align from and to relative to seq (1-based)

		string alnTrace; // descriptive trace info using the characters B, E, M, I, D
	};


	/* static and enum members */
	static const int kNM = 3; // number of matching states
	static const int kNSP = 4; // number of special states
	static const int kNS = kNM + kNSP; // number of total states
	static const string HMM_TAG;
	static const int kMaxProfile = UINT16_MAX + 1;
	static const int kMaxCS = UINT16_MAX + 1;
	static const double kMinGapFrac; // minimum gap fraction comparing to the profile
	static const double CONS_THRESHOLD; // threshold for print upper-case consensus residues
	static const double DEFAULT_ERE; // target mean average relative entropy of the model
	static const Eigen::IOFormat tabFmt;

	/* member functions */
	/* Getters and Setters */
	/**
	 * Get the alphabet
	 */
	const DegenAlphabet* getNuclAbc() const {
		return abc;
	}

	/**
	 * Get the profile name
	 */
	const string& getName() const {
		return name;
	}

	/**
	 * Set the profile name
	 */
	void setName(const string& name) {
		this->name = name;
	}

	/**
	 * Get the profile size
	 */
	int getProfileSize() const {
		return K;
	}

	/**
	 * Set the size of an (un-initialized) profile object
	 * @param size  the designated size of this profile
	 */
	void setProfileSize(int size);

	/**
	 * Set the current profile size to match the K
	 */
	void setProfileSize() {
		return setProfileSize(K);
	}

	/**
	 * Re-calculate the T_MM of B and E states by adding 'wing-retraction'
	 * this method must be called after all T_MM has been read from an hmm file
	 */
	void wingRetract();

	/**
	 * Get the value for the given tag from header options
	 */
	string getOptTag(const string& tag) const;

	/**
	 * Set the header option tag to the given value
	 */
	void setOptTag(const string& tag, const string& val) {
		if(optTags.find(tag) == optTags.end()) // tag not exists
			optTagNames.push_back(tag); // add a new tag
		optTags[tag] = val; // always override
	}

	/**
	 * Set the sequence aligning mode
	 * @param mode  sequence aligning mode, one of GLOBAL, LOCAL, NGCL or CGNL
	 */
	void setSequenceMode(enum align_mode mode);

	/**
	 * Set the special state N and C emission frequencies using given values. No emission for state B and E
	 */
	void setSpEmissionFreq(const Vector4d& freq);

	/**
	 * Set the special state N and C emission frequencies using the embedded background by default.
	 */
	void setSpEmissionFreq() {
		setSpEmissionFreq(hmmBg.getBgEmitPr());
	}

	/** get the CSLen used to build this profile by searching profile2CSIdx */
	int getCSLen() const {
		return L;
	}

	/**
	 * Get the profile location given a index of the original multiple alignment
	 * @param idx  1-based index of the multiple alignment
	 * @return the 1-based position relative to the profile, or 0 is not invalid or unmatched column
	 */
	int getProfileLoc(int idx) const {
		return idx < kMaxProfile ? cs2ProfileIdx[idx] : 0;
	}

	/**
	 * Get the consensus location given a index on the profile
	 * @param idx  1-based index of the BandedHMMP7 profile
	 * @return the 1-based position relative to the consensus sequence, or 0 if out of range
	 */
	int getCSLoc(int idx) const {
		return idx < kMaxCS ? profile2CSIdx[idx] : 0;
	}

	/**
	 * Get the given tag at given loc on the profile
	 * @param tag  optional HMMER3/f tag, i.e. CONS, RF, MM, CS
	 * @param loc  1-based location on the HMM profile
	 * @return  the optional tag value at the match emission lines, or "-" if not defined
	 */
	string getLocOptTag(string tag, int loc) const {
		map<string, vector<string> >::const_iterator it = locOptTags.find(tag);
		if(it == locOptTags.end() || !(loc >= 0 && loc < it->second.size())) /* non-existing tag or out-of-bound loc */
			return "-";
		else
			return it->second[loc];
	}

	/**
	 * Set the given tag at given loc on the profile
	 * @param tag  optional HMMER3/f tag, i.e. CONS, RF, MM, CS
	 * @param val  tag value to set
	 * @param loc  1-based location on the HMM profile
	 */
	void setLocOptTag(string tag, string val, int loc) {
		if(locOptTags[tag].empty())
			locOptTags[tag].resize(K + 1); /* position 0 is dummy */
		locOptTags[tag][loc] = val;
	}

	/**
	 * Prepare the ViterbiScore DP matrices so they are ready for Viterbi algorithm
	 */
	ViterbiScores& prepareViterbiScores(ViterbiScores& vs) const;


	/**
	 * build a known VPath using calculated coordinates
	 * @param csStart  1-based consensus start
	 * @param csEnd  1-based consensus end
	 * @param csFrom  1-based seq start
	 * @param csTo  1-based seq end
	 * @return  a new VPath
	 */
	ViterbiAlignPath buildAlignPath(const CSLoc& csLoc, int csFrom, int csTo) const;

	/**
	 * Calculate full Viterbi DP scores w/o known "seed" alignment region
	 * the S values are updated in the VScore
	 */
	void calcViterbiScores(const PrimarySeq& seq, ViterbiScores& vs) const;

	/**
	 * Calculate banded Viterbi DP scores w/ a given known "seed" alignment region
	 * @return a vector of size L giving the final Viterbi lods of the End state
	 */
	void calcViterbiScores(const PrimarySeq& seq, ViterbiScores& vs, const vector<ViterbiAlignPath>& vpaths) const;

	/**
	 * build the ViterbiTrace matrix, using the B, M, I, D as indicator flags
	 * only the cells in the best score path are filled
	 * @param vs  a ViterbiScores with Viterbi Scores initiated
	 * @param vpath  a ViterbiAlignPath of the DP values by one of the calcViterbiScores methods
	 * @param  a given VTrace to store the result
	 */
	void buildViterbiTrace(const ViterbiScores& vs, ViterbiAlignTrace& vtrace) const;

	/**
	 * Build the global alignment string using calculated scores and backtrace path
 	 * @param seq  the original seq
	 * @param vs  a calculated VScore
	 * @param vtrace  a calcluated VTrace
	 * @return  the global aligned sequence of the query seq, i.e. "AC--GTCGA---ACGNC---";
	 */
	string buildGlobalAlign(const PrimarySeq& seq, const ViterbiScores& vs, const ViterbiAlignTrace& vtrace) const;

	/**
	 * build hmm from a MSA, override any old data
	 * @param msa  a MSA object
	 * @param symfrac  threshold for defining a consensus site in HMM model
	 * @param prior  Dirichlet model based prior models
	 * @param name  optional model name
	 * @return this object
	 */
	BandedHMMP7& build(const MSA& msa, double symfrac,
			const BandedHMMP7Prior& prior, const string& name = "");

	/**
	 * Scale the current model's transition and emission matrix by a constant factor
	 * current model is supposed to be a raw count-based model,
	 * without previous calling of normalize() or estimateParameters()
	 * The new model will also be an invalid model
	 * The invalid model needed to be normalized or estimatedParameters
	 * @param r  constants to scale
	 */
	void scale(double r);

	/**
	 * Normalize an invalid model without using a prior,
	 * the model is usually raw counts or scaled by calling scale()
	 */
	void normalize();

private:
	/* core fields */
	int K; // profile length
	int L; // CSLen used to train this HMM
	/* Transition cost matrices
	 * Note that index 0 indicating B state,
	 * and index K indicating E state
	 */
	vector<Matrix3d> Tmat;

	/* Emission cost matrices */
	Matrix4Xd E_M; /* emission probabilities from Mk node, k = 0, 1, ..., K */
	Matrix4Xd E_I; /* emission probabilities from Ik node, k = 0, 1, ..., K */
	/* No emission from D state */

	Matrix4Xd E_SP; /* Emission probabilities from special_states sp node */
	MatrixXd T_SP; /* log transition probabilities between special states N, B, E, and C */

	/* Entry and exiting probabilities */
	VectorXd entryPr;
	VectorXd exitPr;
	/* By tuning T_SP, entryPr and exitPr probabilities,
	 * we can control the alignment type regarding to the both the profile and sequence.
	 * For 16S rRNA profile-HMM, alignment to profile is always apparently local, as:
	 * entryPr[k] === 1/K
	 * exitPr[k] === 1/K
	 * Alignment with respect to the sequence can be:
	 * global: T_SP(N,N) = T_SP(C,C) = 0
	 * local: T_SP(N,N) = T_SP(C,C) = T_SP(G,G)
	 * or partial global/local as above combinations
	 */

	/* Log transition probabilities, and log emission probabilities, stored as a duplicate copy for speed */
	vector<Matrix3d> Tmat_cost;

	Matrix4Xd E_M_cost;
	Matrix4Xd E_I_cost;

	Matrix4Xd E_SP_cost;
	MatrixXd T_SP_cost;

	VectorXd entryPr_cost;
	VectorXd exitPr_cost;

	/* Banded HMM limits */
	VectorXi gapBeforeLimit; /* Minimum allowed insertions before given position 1..K, with 0 as dummy position */
	VectorXi gapAfterLimit; /* Minimum allowed insertions after given position 1..K, with 0 as dummy position */
	//VectorXi delBeforeLimit; /* Minimum allowed deletions before given position 1..K, with 0 as dummy position */
	//VectorXi delAfterLimit; /* Minimum allowed deletions after given position 1..K, with 0 as dummy position */

	int cs2ProfileIdx[kMaxProfile + 1]; // MAP index from consensus index -> profile index
	int profile2CSIdx[kMaxCS + 1]; // MAP index from profile index -> consensus index

	BandedHMMP7Bg hmmBg; // background HMMP7 profile

	/* information fields */
	string hmmVersion; // version of this hmm file, default is "progName-progVersion"
	string name; // profile name
	const DegenAlphabet* abc; // Nucleotide alphabet
	int nSeq;  // sequence number used in training
	double effN;  // effective sequence number, used with observed counts and Dirichlet prior info in parameter training

	vector<string> optTagNames; // optional tag names in read in order
	map<string, string> optTags; // all HMMER3 optional tag pairs

	map<string, vector<string> > locOptTags; // other profile loc-specific optional tags in the match emission line

	bool wingRetracted;

	/**
	 * Initialize profile transition matrices,
	 * with all prob matrices filled with zero,
	 * and all cost matrices filled with inf
	 */
	void init_transition_params();

	/**
	 * Initialize profile emission cost matrices
	 * with all prob matrices filled with zero,
	 * and all cost matrices filled with inf
	 */
	void init_emission_params();

	/**
	 * Initialize profile special parameters
	 * with all prob matrices filled with zero,
	 * and all cost matrices filled with inf
	 */
	void init_special_params();

	/**
	 * Reset profile transition cost matrices
	 * with all prob matrices filled with zero,
	 * and all cost matrices filled with inf
	 * but keep the profile size K unchanged
	 */
	void reset_transition_params();

	/**
	 * Reset profile emission cost matrices
	 * with all prob matrices filled with zero,
	 * and all cost matrices filled with inf
	 */
	void reset_emission_params();

	/**
	 * Normalize profile transition cost matrices
	 */
//	void normalize_transition_params();

	/**
	 * Normalize profile emission cost matrices
	 * use 0 for normal matrices and -inf for log matrices
	 */
//	void normalize_emission_params();

	/**
	 * Initialize the banded HMM limits as well as their elements
	 */
	void init_limits();

	/**
	 * reset the profile loc index
	 */
	void reset_index();

	/** extend index to maxLen */
	void extend_index();

	/**
	 * set the profile alignment to local mode by setting entry and exit probabilities
	 */
	void enableProfileLocalMode();

	/**
	 * adjust the profile local mode probabilities to accommodate to the learned probabilities
	 */
	void adjustProfileLocalMode();

	/**
	 * Reset all cost matrices by raw probability matrices
	 */
	void resetCostByProb();

	/**
	 * Reset all raw probability matrices by cost matrices
	 */
	void resetProbByCost();

	/**
	 * calculate the mean relative entropy of this model
	 * only match state emission will be used
	 */
	double meanRelativeEntropy() const;

	/**
	 * Re-estimate the parameters using the given prior and current observed frequencies
	 * (usually unnormalzied due to previous call of scale(double)
	 * @param prior  HMM-prior used to estimate
	 */
	void estimateParams(const BandedHMMP7Prior& prior);

	/**
	 * Normalize an invalid model using a Dirichlet model
	 * the model is usually raw counts or scaled by calling scale()
	 */
//	void normalize();


	/* Private utility functions */
	/** Get the minimum of three values */
	static double min(double Vm, double Vi, double Vd) {
		return std::min(Vm, std::min(Vi, Vd));
	}

	/** Get the minimum of four values */
	static double min(double Vb, double Vm, double Vi, double Vd) {
		return std::min(Vb, min(Vm, Vi, Vd));
	}

	/**
	 * Test whether a string is a valid integer format
	 */
	static bool isInteger(const string& s) {
		int i;
		return sscanf(s.c_str(), "%d", &i) == 1;
	}

	/**
	 * trim leading and tailing space of a string
	 */
	static string trim(const string& str, const string& whitespace = " \t");

	/**
	 * decode the p7_state enum to human-readable characters
	 * return null char of not a defined state
	 */
	static char decode(p7_state state);

	/**
	 * encode the human-readable characters to p7_state enum
	 * throw invalid_argument exception if not a valid state
	 *
	 */
	static p7_state encode(char c);

	/**
	 * calculate the distance to the diagnal of a square starting at (from, start)
	 */
	static int diagnalDist(int i, int j, int from, int start) {
		return (i - from) - (j - start);
	}

	/**
	 * Determine the p7 matching state (M, I, D) on a consensus sequence
	 */
	static p7_state determineMatchingState(const int* cs2ProfileIdx, int loc, int8_t base) {
		bool isPos = cs2ProfileIdx[loc] != cs2ProfileIdx[loc - 1];
		return isPos && base >= 0 ? M : isPos && base < 0 ? D : !isPos && base >= 0 ? I : P;
	}

	/* trace back methods to tell which state the current min is coming from */
	/**
	 * four possibility version of whichMin
	 */
	static char whichMin(double probB, double probM, double probI, double probD, const string& states = "BMID") {
		assert(states.length() == 4);
		string::size_type idx = 0;
		double min = inf;
		if(probB < min) {
			idx = 0;
			min = probB;
		}
		if(probM < min) {
			idx = 1;
			min = probM;
		}
		if(probI < min) {
			idx = 2;
			min = probI;
		}
		if(probD < min) {
			idx = 3;
			min = probD;
		}
		/*std::cerr << "probB:" << probB << " probM:" << probM << " probI:" << probI
				<< " probD:" << probD << " min:" << states[idx] << std::endl;*/
		return states[idx];
	}

	/**
	 * three possibility version of whichMin
	 */
	/*static p7_state whichMin(double probM, double probI, double probD, const string& states = "MID") {
		assert(states.length() == 3);
		string::size_type idx = 0;
		double min = inf;
		if(probM > min) {
			idx = 1;
			min = probM;
		}
		if(probI > min) {
			idx = 2;
			min = probI;
		}
		if(probD > min) {
			idx = 3;
			min = probD;
		}
		//std::cerr << "probB:" << probB << " probM:" << probM << " probI:" << probI
		//		<< " probD:" << probD << " min:" << states[idx] << std::endl;
		return encode(states[idx]);
	}*/

	/**
	 * two possibility version of whichMin
	 */
	static char whichMin(double probM, double probID, const string& states) {
		assert(states.length() == 2);
		string::size_type idx = 0;
		double min = inf;
		if(probM < min) {
			idx = 0;
			min = probM;
		}
		if(probID < min) {
			idx = 1;
			min = probID;
		}
		/*std::cerr << "probM:" << probM << " probID:" << probID << " min:" << states[idx] << std::endl;*/
		return states[idx];
	}

	/* convert an hmm coded string to value */
	static double hmmValueOf(const string& s);

	/* print hmm cost values to ostream */
	static ostream& hmmPrintValue(ostream& out, double val);

	static bool yesOrNo2bool(const string& value) {
		return StringUtils::toLower(value) == "yes";
	}

	static string bool2YesOrNo(bool flag) {
		return flag ? "yes" : "no";
	}

	/**
	 * generate padding sequence given the required length and insert (unagiend) sequence
	 * @param L  required length of padding
	 * @param insert  unaligned sequence to pad
	 * @param mode  aligning mode for the insert
	 * @return  a padding sequence trying to use as much of the insert as possible
	 */
	static string getPaddingSeq(int L, const string& insert, char padCh, padding_mode mode);

	/**
	 * generate padding sequence given the required length
	 */
	static string getPaddingSeq(int L, char padCh) {
		return string(L, padCh);
	}

public:
	/** Merge two multiple alignments */
	static string mergeAlign(const string& fwdAln, const string& revAln);

	/** Merge two multiple alignments */
	static string& mergeWith(string& fwdAln, const string& revAln);


//public:
	/* non-member operators */
	/**
	 * utility function for output an alignment path to a human readable string
	 */
	friend ostream& operator<<(ostream& os, const deque<p7_state> path);

	/**
	 * Read a BandedHMMP7 profile from an hmm file
	 */
	friend istream& operator>>(istream& in, BandedHMMP7& hmm);
	/**
	 * Write a BandedHMMP7 profile into a file in hmm format
	 */
	friend ostream& operator<<(ostream& out, const BandedHMMP7& hmm);

	const string& getHmmVersion() const {
		return hmmVersion;
	}

	void setHmmVersion(const string& hmmVersion) {
		this->hmmVersion = hmmVersion;
	}

	friend class RelativeEntropyTargetFunc;

}; /* BandedHMMP7 */

inline std::string BandedHMMP7::getOptTag(const string& tag) const {
	map<string, string>::const_iterator it = optTags.find(tag);
	return it != optTags.end() ? it->second : "";
}

inline void BandedHMMP7::resetCostByProb() {
	/* reset transitions */
	for(int k = 0; k <= K; ++k)
		Tmat_cost[k] = -Tmat[k].array().log();
	/* reset emissions */
	E_M_cost = -E_M.array().log();
	E_I_cost = -E_I.array().log();
}

inline void BandedHMMP7::resetProbByCost() {
	/* reset transitions */
	for(int k = 0; k <= K; ++k)
		Tmat[k] = (-Tmat_cost[k]).array().exp();
	/* reset emissions */
	E_M = (-E_M_cost).array().exp();
	E_I = (-E_I_cost).array().exp();
}

inline char BandedHMMP7::decode(p7_state state) {
	switch(state) {
	case M:
		return 'M';
	case I:
		return 'I';
	case D:
		return 'D';
	case N:
		return 'N';
	case C:
		return 'C';
	case B:
		return 'B';
	case E:
		return 'E';
	case P:
		return P;
	default:
		return CHAR_MAX;
	}
}

inline BandedHMMP7::p7_state BandedHMMP7::encode(char c) {
	switch(c) {
	case 'M':
		return M;
	case 'I':
		return I;
	case 'D':
		return D;
	case 'N':
		return N;
	case 'C':
		return C;
	case 'B':
		return B;
	case 'E':
		return E;
	case 'P':
		return P;
	default:
		throw std::invalid_argument("Invalid state encountered");
	}
}


inline double BandedHMMP7::hmmValueOf(const string& s) {
	return s != "*" ? ::atof(s.c_str()) : inf;
}

inline ostream& hmmPrintValue(ostream& out, double val) {
	return val != EGriceLab::inf ? out << val : out << "*";
}

/**
 * A relative entropy target functor to calculate relative entropy difference
 * between the current status and a given target average information content (in bits)
 */
struct RelativeEntropyTargetFunc : RootFinder::R2RFunc {
	/**
	 * construct a calculator from copies of hmm and prior
	 */
	RelativeEntropyTargetFunc(double ere, const BandedHMMP7& hmm, const BandedHMMP7Prior& prior) :
		ere(ere), hmm(hmm), prior(prior) { }

	/**
	 * virtual destructor, do nothing
	 */
	virtual ~RelativeEntropyTargetFunc() { }

	/**
	 * calculate the relative entropy by scaling the hmm to effN = x
	 * @override  base class abstract method
	 */
	virtual double operator()(double x);

	double ere;
	BandedHMMP7 hmm;
	BandedHMMP7Prior prior;
};

} /* namespace EGriceLab */

#endif /* BANDEDHMMP7_H_ */
