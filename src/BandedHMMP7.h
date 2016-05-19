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
#include <cmath>
#include <cstdio>
#include <stdint.h> /* for fixed size integers */
#include <iostream>
#include "HmmUFOtuConst.h"
#include "BandedHMMP7Bg.h"
#include "BandedHMMCommons.h"
#include "PrimarySeq.h"
#include "MSA.h"

namespace EGriceLab {
using std::string;
using std::istream;
using std::ostream;
using std::vector;
using std::map;
using std::deque;
using Eigen::Matrix3f;
using Eigen::Matrix4Xf;
using Eigen::Matrix4f;

/**
 * Banded plan7 HMM for 16S rRNA profile alignment
 * Similar to HMMER, the global profile includes N, B, M, I, D, E, C states
 * but not the J (joining) state, so repeated and multiple mapping of same region not allowed
 */
class BandedHMMP7 {
public:
	/* constructors */
	/**
	 * Construct a NULL BandedHMMP7 profile with name as "NULL",
	 * length 0, and values not initialized
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
	enum p7_state { M, I, D, N, C, B, E };

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
		MatrixXf DP_M;  /* the log-odds scores of the best path matching the subsequence X1..i
							to the profile submodel up to the column j, ending with xi being emitted by Mj*/
		MatrixXf DP_I;  /* the log-odds scores of the best path matching the subsequence Xi..i
							to the profile submodel up to the ending in xi being emitted by Ij */
		MatrixXf DP_D;  /* the log-odds score of the best path ending in Dj, and xi being the last character emitted before Dj).*/

		MatrixXf S;     /*  the Dynamic-Programming scoring matrix storing the whole ViterbiScores (log-odds)
		of all combinations of i in 0..L and j in 0..K */

		/* Forward algorithm probability matrices */
		/*MatrixXf F_M;*/  /* forward probability of all alignments up to xi and ends in Mj */
		/*MatrixXf F_I;*/  /* forward probability of all alignments up to xi and ends in Ij */
		/*MatrixXf F_D;*/ /* forward probability of all alignments up to xi and ends in Dj */

		/* Backward algorithm probability matrices*/
		/*MatrixXf B_M;*/ /* backward probability of all alignments up to xi and ends in Mj */
		/*MatrixXf B_I;*/ /* backward probability of all alignments up to xi and ends in Ij */
		/*MatrixXf B_D;*/ /* backward probability of all alignments up to xi and ends in Dj */

		/* Forward-backward scaling vectors to prevent numeric underflow */
/*		VectorXf SV_M;
		VectorXf SV_I;
		VectorXf SV_D;*/

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
	static const float inf; // inf
	static const float infV; // -inf
	static const int kNM; // number of matching states
	static const int kNSP; // number of special states
	static const int kNS; // number of total states
	static const string HMM_TAG;
	static const int kMaxProfile = UINT16_MAX + 1;
	static const int kMaxCS = UINT16_MAX + 1;
	static const float kMaxGapFrac; // maximum gap fraction comparing to the profile

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
	string getHeaderOpt(const string& tag) const {
		map<string, string>::const_iterator it = headOpts.find(tag);
		if(it != headOpts.end()) // tag exists
			return it->second;
		else
			return ""; // use empty string by default
	}

	/**
	 * Set the header option tag to the given value
	 */
	void setHeaderOpt(const string& tag, const string& val) {
		if(headOpts.find(tag) == headOpts.end()) // tag not exists
			headNames.push_back(tag); // add a new tag
		headOpts[tag] = val;
	}

	/**
	 * Set the sequence aligning mode
	 * @param mode  sequence aligning mode, one of GLOBAL, LOCAL, NGCL or CGNL
	 */
	void setSequenceMode(enum align_mode mode);

	/**
	 * Set the special state N and C emission frequencies. No emission for state B and E
	 */
	void setSpEmissionFreq(const Vector4f& freq);

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
	 * @param tag  ptional HMMER3/f tag, i.e. CONS, RF, MM, CS
	 * @param loc  1-based location on the HMM profile
	 * @return  the optional tag value at the match emission lines, or "-" if not defined
	 */
	string getLocOptTag(string tag, int loc) const {
		map<string, map<int, string> >::const_iterator it = locOptTags.find(tag);
		if(it == locOptTags.end() || it->second.find(loc) == it->second.end()) /* non-existing tag or out-of-bound loc */
			return "-";
		else
			return it->second.find(loc)->second;
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
	 * @param csEnd  1-based end on profile on consensus of MSA
	 * @param csFrom  1-based start on sequence
	 * @param csTo  1-based end on sequence
	 * @param cs the matching sequence, including gaps
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
	float buildViterbiTrace(const ViterbiScores& vs, ViterbiAlignPath& vpath) const;

	/**
	 * Build the global alignment string using calculated scores and backtrace path
	 * @param vs  a ViterbiScores with Viterbi Scores initiated
	 * @param vpath  a ViterbiAlignPath of the DP values by one of the calcViterbiScores methods
	 * @return  the global aligned sequence of the query seq, i.e. "AC--GTCGA---ACGNC---";
	 */
	string buildGlobalAlignSeq(const ViterbiScores& vs, const ViterbiAlignPath& vpath) const;

	/* static member methods */
	static BandedHMMP7* build(const MSA* msa, double symfrac, const string& name = "unnamed");

private:
	string version; // version of the program generated this hmm file, default is progName-progVersion
	string name; // profile name
	const DegenAlphabet* abc; // Nucleotide alphabet
	int K; // profile length
	BandedHMMP7Bg hmmBg; // background HMMP7 profile
	vector<string> headNames; // all header names, if set
	map<string, string> headOpts; // header options the profile, as compatitable to HMMER3/b
	int cs2ProfileIdx[kMaxProfile + 1]; // MAP index from consensus index -> profile index
	int profile2CSIdx[kMaxCS + 1]; // MAP index from profile index -> consensus index
	map<string, map<int, string> > locOptTags; // other profile loc-specific optional tags in the match emission line

	/* Transition probability matrices */
	/*
	 * Note that index 0 indicating B state,
	 * and index 1 indicating E state
	 */
	vector<Matrix3f> Tmat;

	/* Emission probability matrices */
	Matrix4Xf E_M; /* emission probabilities from Mk node */
	Matrix4Xf E_I; /* emission probabilities from Ik node */
	/* No emission from D state */

	Matrix4Xf E_SP; /* Emission probabilities from special_states sp node */
	MatrixXf T_SP; /* log transition probabilities between special states N, B, E, and C */

	/* Entry and exiting probabilities */
	VectorXf entryPr;
	VectorXf exitPr;
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
	vector<Matrix3f> Tmat_log;

	Matrix4Xf E_M_log;
	Matrix4Xf E_I_log;

	Matrix4Xf E_SP_log;
	MatrixXf T_SP_log;

	VectorXf entryPr_log;
	VectorXf exitPr_log;

	/* Special wing retracted T_MM transition probabilities,
	 * by adding the B->D1->...->Dk-1 to the B->Mk probability
	 * and adding the Dk+1->Dk+2->...->E to the Mk->E probabitity
	 */
	//MatrixXf T_MM_retract;
	//MatrixXf T_MM_retract_log;

	/* Banded HMM limits */
	VectorXi gapBeforeLimit; /* Maximum allowed insertions before given position 1..K, with 0 as dummy position */
	VectorXi gapAfterLimit; /* Maximum allowed insertions after given position 1..K, with 0 as dummy position */
	//VectorXi delBeforeLimit; /* Maximum allowed deletions before given position 1..K, with 0 as dummy position */
	//VectorXi delAfterLimit; /* Maximum allowed deletions after given position 1..K, with 0 as dummy position */

	/* Initialization flags */
	bool transInit;
	bool emisInit;
	bool spInit;
	bool limitInit;
	bool indexInit;
	bool wingRetracted;

	/**
	 * Initialize profile transition probability matrices, but not their elements
	 */
	void init_transition_params();

	/**
	 * Initialize profile emission probability matrices, but not their elements
	 */
	void init_emission_params();

	/**
	 * Initialize profile special parameters, but not their elements
	 */
	void init_special_params();

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
	 * test whether given viterbi align path is valid
	 * @param L  query length
	 * @param start  1-based index on profile
	 * @param end  1-based index on profile
	 * @param from 1-based index on query
	 * @param to 1-based index on query
	 * @param path  path string
	 */
	bool isValidAlignPath(const ViterbiAlignPath& vpath) const;

	/* Private static utility functions */
	/** Get the maximum of three values */
	static float max(float Vm, float Vi, float Vj) {
		return std::max(Vm, std::max(Vi, Vj));
	}

	/** Get the maximum of four values */
	static float max(float Vb, float Vm, float Vi, float Vj) {
		return std::max(Vb, max(Vm, Vi, Vj));
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

	/* trace back methods to tell which state the current max is coming from */
	/**
	 * four possibility version of whichMax
	 */
	static char whichMax(float probB, float probM, float probI, float probD, const string& states = "BMID") {
		assert(states.length() == 4);
		string::size_type idx = 0;
		float max = infV;
		if(probB > max) {
			idx = 0;
			max = probB;
		}
		if(probM > max) {
			idx = 1;
			max = probM;
		}
		if(probI > max) {
			idx = 2;
			max = probI;
		}
		if(probD > max) {
			idx = 3;
			max = probD;
		}
		/*std::cerr << "probB:" << probB << " probM:" << probM << " probI:" << probI
				<< " probD:" << probD << " max:" << states[idx] << std::endl;*/
		return states[idx];
	}

	/**
	 * three possibility version of whichMax
	 */
	/*static p7_state whichMax(float probM, float probI, float probD, const string& states = "MID") {
		assert(states.length() == 3);
		string::size_type idx = 0;
		float max = infV;
		if(probM > max) {
			idx = 1;
			max = probM;
		}
		if(probI > max) {
			idx = 2;
			max = probI;
		}
		if(probD > max) {
			idx = 3;
			max = probD;
		}
		//std::cerr << "probB:" << probB << " probM:" << probM << " probI:" << probI
		//		<< " probD:" << probD << " max:" << states[idx] << std::endl;
		return encode(states[idx]);
	}*/

	/**
	 * two possibility version of whichMax
	 */
	static char whichMax(float probM, float probID, const string& states) {
		assert(states.length() == 2);
		string::size_type idx = 0;
		float max = -std::numeric_limits<float>::infinity();
		if(probM > max) {
			idx = 0;
			max = probM;
		}
		if(probID > max) {
			idx = 1;
			max = probID;
		}
		/*std::cerr << "probM:" << probM << " probID:" << probID << " max:" << states[idx] << std::endl;*/
		return states[idx];
	}

	/* convert an hmm coded string to value */
	static float hmmValueOf(const string& s);

	/* non-member operators */
	/**
	 * utility function for output an alignment path to a human readable string
	 */
	friend ostream& operator<<(ostream& os, const deque<p7_state> path);

	/**
	 * Read a BandedHMMP7 profile from an hmm file
	 */
	friend istream& operator>>(istream& is, BandedHMMP7& hmm);
	/**
	 * Write a BandedHMMP7 profile into a file in hmm format
	 */
	friend ostream& operator<<(ostream& os, const BandedHMMP7& hmm);
}; /* BandedHMMP7 */

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
	default:
		return '\0';
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
	default:
		throw std::invalid_argument("Invalid state encountered");;
	}
}

inline float BandedHMMP7::hmmValueOf(const string& s) {
	return s != "*" ? ::atof(s.c_str()) : BandedHMMP7::inf;
}

} /* namespace EGriceLab */

#endif /* BANDEDHMMP7_H_ */
