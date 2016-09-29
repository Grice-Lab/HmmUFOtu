/*
 * MSA.h
 *
 *  Created on: Jul 23, 2015
 *      Author: zhengqi
 */

#ifndef MSA_H_
#define MSA_H_
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include "SeqCommons.h"
#include "StringUtils.h"
#include "PrimarySeq.h"
#include "DigitalSeq.h"

namespace EGriceLab {
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::out_of_range;
using std::invalid_argument;
using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::Map;

/**
 * A class for Multiple Sequence Alignment
 * All aligned sequences are stored in concatenated form to save space
 * @version v1.1
 * @since v1.1
 */
class MSA {
public:
	/* constructors */
	/** destructor, do nothing */
	virtual ~MSA() {
		clear();
	}

	/** Setters and Getters */
	const string& getAlphabet() const {
		return alphabet;
	}

	const DegenAlphabet* getAbc() const {
		return abc;
	}

	const string& getName() const {
		return name;
	}

	void setName(const string& name) {
		this->name = name;
	}

	const string& getConcatMsa() const {
		return concatMSA;
	}

	const vector<string>& getSeqNames() const {
		return seqNames;
	}

	/*
	 * Get the consensus seq of this MSA
	 * @return the consensus seq
	 */
	const string& getCS() const {
		return CS;
	}

	unsigned getCSLen() const {
		return csLen;
	}

	unsigned getNumSeq() const {
		return numSeq;
	}

	unsigned long getTotalNumGap() const {
		return gapCount.sum();
	}

	/* member methods */
	/**
	 * Prune this MSA by removing all gapped site
	 * @return the modified MSA object
	 */
	MSA& prune();

	/**
	 * get the total length of this MSA
	 * @return the total MSA length w/ gaps
	 */
	unsigned long getMSALen() const;

	/**
	 * get the total length of non-gapped residuals of this MSA
	 * @return the total MSA length wo/ gaps
	 */
	unsigned long getMSANonGapLen() const;

	/**
	 * get the ith seq name
	 * @param i  seq position
	 * @return seq name of the ith seq
	 */
	const string& seqNameAt(unsigned i) const;

	/**
	 * Get the residual at given seq and CS pos
	 * @param i  seq position
	 * @param j  CS position
	 * @return  residual at this pos
	 */
	char residualAt(unsigned i, unsigned j) const;

	/**
	 * Get the encoded residual at given seq and CS pos
	 * @param i  seq position
	 * @param j  CS position
	 * @return  encoded value between 0..alphabetSize - 1 at this pos
	 */
	int8_t encodeAt(unsigned i, unsigned j) const;

	/**
	 * Test whether a given residual is a gap
	 * @param i  seq position
	 * @param j  CS position
	 * @return  gap or not at this residual
	 */
	bool isGapAt(unsigned i, unsigned j) const {
		return abc->isGap(residualAt(i, j));
	}

	/**
	 * Get the given aligned seq
	 * @param i  seq position
	 * @return the i-th seq
	 */
	string seqAt(unsigned i) const;

	/**
	 * Get the PrimarySeq at given pos
	 * @param i  seq position
	 * @return the i-th PrimarySeq
	 */
	PrimarySeq primarySeqAt(unsigned i) const;

	/**
	 * Get the DigitalSeq at given pos
	 * @param i  seq position
	 * @return the i-th DigitalSeq
	 */
	DigitalSeq dsAt(unsigned i) const;

	/**
	 * Get the alignment string at given pos
	 * @param j  CS position
	 * @return  the alignment at the j-th CS pos
	 */
	string alignAt(unsigned j) const;

	/**
	 * Get the seq start on CS
	 * @param i  seq index
	 * @return  ith seq start
	 */
	int seqStart(unsigned i) const;

	/**
	 * Get the seq end on CS
	 * @param i  seq index
	 * @return  ith seq start
	 */
	int seqEnd(unsigned i) const;

	/**
	 * Get the seq length on CS
	 * @param i  seq index
	 * @return  ith seq ungapped length
	 */
	int seqLength(unsigned i) const;

	/**
	 * Get the consensus residual at given pos
	 * @param j  CS position
	 * @return the consensus residual as the most frequent one
	 */
	char CSResidualAt(unsigned j) const;

	/**
	 * Get the consensus non-gap residual at given pos
	 * @param j  CS position
	 * @return the consensus base as the most frequent one
	 */
	char CSBaseAt(unsigned j) const;

	/**
	 * Get the identity of aligned residuals at given pos
	 * @param j  CS position
	 * @return  identity that is the frection of residuals that match to the consensus one
	 */
	double identityAt(unsigned j) const;

	/**
	 * Get the weighted identity of aligned residuals at given pos
	 * @param j  CS position
	 * @return  weighted identity that is the fraction of residuals that match to the consensus one
	 */
	double wIdentityAt(unsigned j) const;

	/**
	 * Get the weight of the ith seq
	 * @param i  seq index
	 * @return  weight of the ith seq
	 */
	double getSeqWeight(unsigned i) const;

	/**
	 * Get the fraction of gaps of at given pos
	 * @param j  CS position
	 * @return  fraction of gaps at this pos
	 */
	double gapFrac(unsigned j) const;

	/**
	 * Get the fraction of symbols of at given pos
	 * @param j  CS position
	 * @return  fraction of symbols at this pos
	 */
	double symFrac(unsigned j) const;

	/**
	 * Get the fraction of symbols of at given pos
	 * @param j  CS position
	 * @return  weighted fraction of symbols at this pos
	 */
	double symWFrac(unsigned j) const;

	/**
	 * Get the symbol frequency at given pos
	 * @param j  CS position
	 * @return  frequency count of symbols at this pos
	 */
	VectorXd symFreq(unsigned j) const;

	/**
	 * Get the weighted symbol frequency at given pos
	 * @param j  CS position
	 * @return  weighted frequency count of symbols at this pos
	 */
	VectorXd symWFreq(unsigned j) const;

	/**
	 * Get the relative entropy (or information content) at given pos
	 * @param j  CS position
	 * @return  entropy calculated with weighted count
	 */
	double entropyAt(unsigned j) const;

	/**
	 * Get the residual frequency of each base of the entire MSA
	 */
	VectorXd resFreq() const;

	/**
	 * Get the weighted residual frequency of each base of the entire MSA
	 */
	VectorXd resWFreq() const;

	/**
	 * Get the effective number of sequence of this MSA as the sum of the seq weights
	 */
	double getEffectSeqNum() const;

	/**
	 * Update the count matrices of this object, also update auxiliary indices
	 */
	void updateRawCounts();

	/**
	 * Update entropy of the PSSM
	 */
	void updateEntropy();

	/**
	 * Update the seqWeight of this object
	 */
	void updateSeqWeight();

	/**
	 * Update the count matrices of this object
	 */
	void updateWeightedCounts();

	/**
	 * Normalize seqWeights so the the average bit per consensus site is set to DEFAULT_ENTROPY_PER_BASE
	 * @return  the normalization factor z
	 */
	double normalizeSeqWeight(double ere = DEFAULT_ENTROPY_PER_BASE, double symfrac = DEFAULT_CONSENSUS_FRAC);

	/**
	 * Save this MSA object to a binary file
	 * @param f  the binary output
	 * @return  true if everything is good after saved
	 */
	std::ostream& save(std::ostream& out);

	/**
	 * Save this MSA object to a given file
	 * @param filename  filename to save to
	 * @param format  MSA file format
	 */
	bool saveMSAFile(const string& filename, const string& format);

	/**
	 * Save this MSA object to a file in FASTA format
	 * @param filename  filename to save to
	 */
	bool saveFastaFile(const string& filename);


	/* static methods */
	/**
	 * Load an MSA file in given format
	 * @param filename  MSA file name
	 * @param format  MSA file format
	 * @return  a newly constructed MSA pointer
	 * @throws invalid_argument if alphabet or format is not known
	 */
	static MSA* load(std::istream& in);

	/**
	 * Load an MSA binary file
	 * @param f  the binary input
	 * @return  a newly constructed MSA pointer
	 * @throws invalid_argument if alphabet or format is not known
	 */
	static MSA* loadMSAFile(const string& alphabet, const string& filename, const string& format);

	/**
	 * Load an MSA file in fasta format
	 * @param filename  MSA file name
	 * @return  a newly constructed MSA pointer
	 * @throws invalid_argument if alphabet or format is not known
	 */
	static MSA* loadFastaFile(const string& alphabet, const string& filename);

	/**
	 * calculate the relative entropy in bits using given observed frequency and the background frequency
	 */
	static double relEntropy(const VectorXd& p, const VectorXd& bg);

	/**
	 * calculate the absolute entropy in bits using given observed frequency
	 */
	static double absEntropy(const VectorXd& p);

private:
	/* constructors */
	/**
	 * construct an MSA with given alphabet
	 * @throw invalid_argument if the alphabet is not known
	 */
	explicit MSA(const string& alphabet = "dna") : alphabet(alphabet), abc(SeqCommons::getAlphabetByName(alphabet)),
		numSeq(0), csLen(0), isPruned(false),
		startIdx(NULL), lenIdx(NULL), endIdx(NULL),
		resCountBuf(NULL), gapCountBuf(NULL),
		seqWeightBuf(NULL), posEntropyBuf(NULL),
		resWCountBuf(NULL), gapWCountBuf(NULL),
		resCount(NULL, 0, 0), gapCount(NULL, 0),
		seqWeight(NULL, 0), posEntropy(NULL, 0),
		resWCount(NULL, 0, 0), gapWCount(NULL, 0)
	{  }

	/* Disable copy and assignment constructors */
	MSA(const MSA& other);
	MSA& operator=(const MSA& other);

	/* Clear the heap memories */
	void clear();

	/* reset count values */
	void resetRawCount();

	/* reset entropy */
	void resetEntropy();

	/* reset seq weights */
	void resetSeqWeight();

	/* reset weighted count values */
	void resetWeightedCount();

	/* calculate the CS if not provided by the MSA file */
	void calculateCS();

	string alphabet;
	const DegenAlphabet* abc;
	string name;
	unsigned numSeq; /* number of sequences */
	unsigned csLen;  /* consensus seq length */
	vector<string> seqNames; /* seq names stored in their occurring order */
	//vector<string> ids;
	string concatMSA; // concatenated MSA
	string CS;        // Consensus Sequence
	bool isPruned; // flag for whether this MS is pruned
	/* auxiliary data to remember each sequence start, end and length (non-gapped) */
	int* startIdx; /* start position on CS */
	int* endIdx; /* end position on CS */
	int* lenIdx; /* unmapped length */

	/* raw matrix/vector buffers for residual & gap count */
	int* resCountBuf; /* raw buffer for Residual count */
	int* gapCountBuf; /* raw buffer for gap count */
	double* seqWeightBuf; /* raw buffer for sequence weight */
	double* posEntropyBuf; /* raw buffer for position entropy (in bits) */
	double* resWCountBuf; /* raw buffer for weighted residual count */
	double* gapWCountBuf; /* raw buffer for weighted gap count */

	/* matrix/vector for residual & gap count */
	Map<MatrixXi> resCount; /* Residual count matrix w/ alphabet-size X CSLen dimension */
	Map<VectorXi> gapCount; /* gap count array w/ CSLen length */
	Map<VectorXd> seqWeight; /* Sequence weight for each seq w/ numSeq length */
	Map<VectorXd> posEntropy; /* Position entropy of the PSSM (in bits) */
	Map<MatrixXd> resWCount; /* weighted residual count matrix w/ alphabet-size X CSLen dimension */
	Map<VectorXd> gapWCount; /* weighted gap count array w/ CSLen length */

	/* static members */
public:
	static const double DEFAULT_ENTROPY_PER_BASE = 1;
	static const double NAT2BIT;
	static const double DEFAULT_CONSENSUS_FRAC = 0.5;
};

inline unsigned long MSA::getMSALen() const {
	return static_cast<unsigned long> (numSeq) * static_cast<unsigned long> (csLen);
}

inline unsigned long MSA::getMSANonGapLen() const {
	return getMSALen() - getTotalNumGap();
}

inline const string& MSA::seqNameAt(unsigned i) const {
	return seqNames[i];
}

inline char MSA::residualAt(unsigned i, unsigned j) const {
/*	if(!(i >= 0 && i < numSeq && j >= 0 && j < csLen))
		throw out_of_range("residual is out of range");*/
	return concatMSA[i * csLen + j];
}

inline int8_t MSA::encodeAt(unsigned i, unsigned j) const {
	return abc->encode(::toupper(concatMSA[i * csLen + j]));
}

inline string MSA::seqAt(unsigned i) const {
	return concatMSA.substr(i * csLen, csLen);
}

inline PrimarySeq MSA::primarySeqAt(unsigned i) const {
	return PrimarySeq(abc, seqNameAt(i), seqAt(i));
}

inline DigitalSeq MSA::dsAt(unsigned i) const {
	return DigitalSeq(abc, seqNameAt(i), seqAt(i));
}

inline string MSA::alignAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // only check j range once
		throw out_of_range("CS pos is out of range");
	string aln;
	for(unsigned i = 0; i < numSeq; ++i)
		aln.push_back(concatMSA[i * csLen + j]);
	return aln;
}

inline int MSA::seqStart(unsigned i) const {
	return startIdx[i];
}

inline int MSA::seqEnd(unsigned i) const {
	return endIdx[i];
}

inline int MSA::seqLength(unsigned i) const {
	return lenIdx[i];
}

inline double MSA::getSeqWeight(unsigned i) const {
	return seqWeight(i);
}

inline VectorXd MSA::symFreq(unsigned j) const {
	return resCount.col(j).cast<double>();
}

inline VectorXd MSA::symWFreq(unsigned j) const {
	return resWCount.col(j);
}

inline double MSA::entropyAt(unsigned j) const {
	return posEntropy(j);
}

inline VectorXd MSA::resFreq() const {
	VectorXd freq = resCount.rowwise().sum().cast<double>();
	return freq / freq.sum();
}

inline VectorXd MSA::resWFreq() const {
	VectorXd freq = resWCount.rowwise().sum();
	return freq / freq.sum();
}

inline double MSA::getEffectSeqNum() const {
	return seqWeight.sum();
}

inline MSA* MSA::loadMSAFile(const string& alphabet,
		const string& filename, const string& format) {
	MSA* msa = NULL;
	if(format == "fasta")
		msa = loadFastaFile(alphabet, filename);
	else
		throw invalid_argument("Unsupported MSA file format '" + format + "'");
	/* construct the CS, if it is not set by the MSA file yet */
	msa->calculateCS();
	return msa;
}

inline bool MSA::saveMSAFile(const string& filename, const string& format) {
	if(format == "fasta")
		return MSA::saveFastaFile(filename);
	else throw invalid_argument("Cannot save MSA to file, unsupported MSA file format " + format);
}

inline double MSA::relEntropy(const VectorXd& p, const VectorXd& bg) {
	assert(p.rows() == bg.rows());
	double entropy = 0;
	for(VectorXd::Index i = 0; i < p.rows(); ++i)
		if(p(i) > 0)
			entropy += p(i) * ::log(static_cast<double> (p(i)) / static_cast<double> (bg(i)));
	return entropy * NAT2BIT;
}

inline double MSA::absEntropy(const VectorXd& p) {
	double entropy = 0;
	for(VectorXd::Index i = 0; i < p.rows(); ++i)
		if(p(i) > 0)
			entropy += p(i) * ::log(static_cast<double> (p(i)));
	return entropy * NAT2BIT;
}

} /* namespace EGriceLab */

#endif /* MSA_H_ */
