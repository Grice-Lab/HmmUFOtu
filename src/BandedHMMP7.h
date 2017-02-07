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
#include <eigen3/Eigen/Dense>
#include <limits>
#include <cmath>
#include <cstdio>
#include <climits>
#include <stdint.h> /* for fixed size integers */
#include <iostream>
#include "HmmUFOtuConst.h"
#include "BandedHMMP7Bg.h"
#include "BandedHMMP7Prior.h"
#include "BandedHMMCommons.h"
#include "StringUtils.h"
#include "PrimarySeq.h"
#include "DigitalSeq.h"
#include "MSA.h"
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
	 * Construct a BandedHMMP7 with given length
	 */
	BandedHMMP7(const string& name, int K, const DegenAlphabet* abc = SeqCommons::nuclAbc);

	/* nested members */
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

	enum align_mode {
		GLOBAL,
		LOCAL,
		NGCL /* N' global C' local */,
		CGNL /* C' global N' local */
	};

	/**
	 * A nested class storing public accessible Viterbi Scoring matrices
	 */
	class ViterbiScores {
	public:
		/* constructors*/
		/* default constructor to construct a ViterbiScores object without any given seq, do not do anything with it */
		ViterbiScores() : seq(NULL) { }

		/* construct a ViterbiScores object with given DigitalSeq */
		explicit ViterbiScores(const PrimarySeq& seq) : seq(&seq) { }

		const PrimarySeq* seq; // a pointer to a const DigitalSeq
		/* Viterbi log-odd scoring matrices */
		MatrixXd DP_M;  /* the cost of the best path matching the subsequence X1..i
							to the profile submodel up to the column j, ending with xi being emitted by Mj*/
		MatrixXd DP_I;  /* the cost of the best path matching the subsequence Xi..i
							to the profile submodel up to the ending in xi being emitted by Ij */
		MatrixXd DP_D;  /* the lcost of the best path ending in Dj, and xi being the last character emitted before Dj).*/

		MatrixXd S;     /*  the Dynamic-Programming scoring matrix storing the whole ViterbiScores (log-odds)
		of all combinations of i in 0..L and j in 0..K */

		/* Forward algorithm cost matrices */
		/*MatrixXd F_M;*/  /* forward cost of all alignments up to xi and ends in Mj */
		/*MatrixXd F_I;*/  /* forward cost of all alignments up to xi and ends in Ij */
		/*MatrixXd F_D;*/ /* forward cost of all alignments up to xi and ends in Dj */

		/* Backward algorithm cost matrices*/
		/*MatrixXd B_M;*/ /* backward cost of all alignments up to xi and ends in Mj */
		/*MatrixXd B_I;*/ /* backward cost of all alignments up to xi and ends in Ij */
		/*MatrixXd B_D;*/ /* backward cost of all alignments up to xi and ends in Dj */

		/* Forward-backward scaling vectors to prevent numeric underflow */
/*		VectorXd SV_M;
		VectorXd SV_I;
		VectorXd SV_D;*/

	};

	class ViterbiAlignPath {
	public:
		int L; // sequence length
		int K; // profile length
		int N; // # of known path segments
		vector<int> start, end; // known align path coordinates relative to profile (1-based)
		vector<int> from, to; // known align path relative to sequence (1-based)
		vector<int> numIns, numDel; // number of insertion and deletions in this path
		int alnStart, alnEnd; // final align start and end (1-based)
		int alnFrom, alnTo; // final align from and to relative to seq (1-based)
		//deque<p7_state> path;
		string alnPath; // final align path
		/** Back-trace matrix indicating which states this state is from
		 * -1 indicates not determined yet
		 */
		//MatrixXi TRACE;
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
	 * Set the special state N and C emission frequencies. No emission for state B and E
	 */
	void setSpEmissionFreq(const Vector4d& freq);

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
	 * Initialize the HMM algorithm scoring schema without a given sequence,
	 * nothing is initialized, has same effect as the default constructor,
	 * but kept for format consistency
	 */
	ViterbiScores initViterbiScores() const {
		return ViterbiScores();
	}

	/**
	 * Initialize the HMM algorithm scoring schema with a given sequence,
	 * set the all members and matrix elements to default values
	 */
	ViterbiScores initViterbiScores(const PrimarySeq& seq) const;

	/**
	 * Reset the ViterbiScore matrices size and elements to fit the new seq
	 */
	ViterbiScores& resetViterbiScores(ViterbiScores& vs, const PrimarySeq& seq) const;

	/**
	 * Initialize a Viterbi Align Path for the HMM algorithm
	 * set all members and matrix elements to default values
	 */
	ViterbiAlignPath initViterbiAlignPath(int L) const;

	/**
	 * Set the Viterbi Align Path to fit the new length
	 */
	ViterbiAlignPath& resetViterbiAlignPath(ViterbiAlignPath& vpath, int L) const;

	/**
	 * Set the forward scores in ViterbiScores to default values
	 */
/*	void resetForwardBackwardScores(ViterbiScores& bHmm) const;*/

	/**
	 * Set the tracke scores in ViterbiScores to default values
	 */
/*	void resetTraceScores(ViterbiScores& bHmm) const;*/

	/**
	 * Set a known alignment path (from database searches) for the BandedHMMP7 dynamic programming
	 * fill the TRACE matrix with states
	 * @param vpath  a ViterbiAlignPath
	 * @param csStart  1-based start on consensus of MSA
	 * @param csEnd  1-based end on consensus of MSA
	 * @param csFrom  1-based start on sequence
	 * @param csTo  1-based end on sequence
	 */
	void addKnownAlignPath(ViterbiAlignPath& vpath, const CSLoc& csLoc, int csFrom, int csTo) const;

	/**
	 * Calculate full Viterbi DP scores w/o known "seed" alignment region
	 * @return a vector of size L giving the final Viterbi lods of the End state
	 */
	void calcViterbiScores(ViterbiScores& vs) const;

	/**
	 * Calculate banded Viterbi DP scores w/ a given known "seed" alignment region
	 * @return a vector of size L giving the final Viterbi lods of the End state
	 */
	void calcViterbiScores(ViterbiScores& vs, const ViterbiAlignPath& vpath) const;

	/**
	 * build the ViterbiTrace matrix, using the B, M, I, D as indicator flags
	 * only the cells in the best score path are filled
	 * @param vs  a ViterbiScores with Viterbi Scores initiated
	 * @param vpath  a ViterbiAlignPath of the DP values by one of the calcViterbiScores methods
	 * @return  a track string of the entire sequence match, with same length as profile length
	 */
	double buildViterbiTrace(const ViterbiScores& vs, ViterbiAlignPath& vpath) const;

	/**
	 * Build the global alignment string using calculated scores and backtrace path
	 * @param vs  a ViterbiScores with Viterbi Scores initiated
	 * @param vpath  a ViterbiAlignPath of the DP values by one of the calcViterbiScores methods
	 * @return  the global aligned sequence of the query seq, i.e. "AC--GTCGA---ACGNC---";
	 */
	string buildGlobalAlignSeq(const ViterbiScores& vs, const ViterbiAlignPath& vpath) const;

	/**
	 * Build the global aligned DigitalSeq using calculated scores and backtrace path
	 * @param vs  a ViterbiScores with Viterbi Scores initiated
	 * @param vpath  a ViterbiAlignPath of the DP values by one of the calcViterbiScores methods
	 * @return  the global aligned sequence of the query seq
	 */
	DigitalSeq buildGlobalAlignDS(const ViterbiScores& vs, const ViterbiAlignPath& vpath) const;

	/* static member methods */
	static BandedHMMP7 build(const MSA& msa, double symfrac,
			const BandedHMMP7Prior& prior, const string& name = "unnamed");

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

	/* Initialization flags */
	bool transInit;
	bool emisInit;
	bool spInit;
	bool limitInit;
	bool indexInit;
	bool wingRetracted;

	/**
	 * Initialize profile transition cost matrices, but not their elements
	 */
	void init_transition_params();

	/**
	 * Initialize profile emission cost matrices, but not their elements
	 */
	void init_emission_params();

	/**
	 * Initialize profile special parameters, but not their elements
	 */
	void init_special_params();

	/**
	 * Reset profile transition cost matrices
	 * use 0 for normal matrices and -inf for log matrices
	 */
	void reset_transition_params();

	/**
	 * Reset profile emission cost matrices
	 * use 0 for normal matrices and -inf for log matrices
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
	 * Initialize the profile loc index
	 */
	void init_index();

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

	/**
	 * test whether given viterbi align path is valid
	 * @param L  query length
	 * @param start  1-based index on profile
	 * @param end  1-based index on profile
	 * @param from 1-based index on query
	 * @param to 1-based index on query
	 * @param path  path string
	 */
	bool isValidAlignPath(const ViterbiAlignPath& vpath) const;

	/* Private utility functions */
	/** Get the minimum of three values */
	static double min(double Vm, double Vi, double Vj) {
		return std::min(Vm, std::min(Vi, Vj));
	}

	/** Get the minimum of four values */
	static double min(double Vb, double Vm, double Vi, double Vj) {
		return std::min(Vb, min(Vm, Vi, Vj));
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

/*	*
	 * Determine the p7 special state (M, I, D) on a consensus sequence

	static p7_state determineSpecialState(int csStart, int csEnd, int loc) {
		return loc < csStart ? B : loc > csEnd ? E : M;
	}*/

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
		throw std::invalid_argument("Invalid state encountered");;
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
