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
#include "SeqCommons.h"

namespace EGriceLab {
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::out_of_range;
using std::invalid_argument;

/**
 * A class for Multiple Sequence Alignment
 * All aligned sequences are stored in concatenated form to save space
 * @version v1.1
 * @since v1.1
 */
class MSA {
public:
	/** Setters and Getters */
	const string& getAlphabet() const {
		return alphabet;
	}

	const DegenAlphabet* getAbc() const {
		return abc;
	}

	const string& getConcatMsa() const {
		return concatMSA;
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
		return totalNumGap;
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
	 * Get the residual at given seq and CS pos
	 * @param i  seq position
	 * @param j  CS position
	 * @return  residual at this pos
	 */
	char residualAt(unsigned i, unsigned j) const;

	/**
	 * Get the decoded residual at given seq and CS pos
	 * @param i  seq position
	 * @param j  CS position
	 * @return  decoded value between 0..alphabetSize - 1 at this pos
	 */
	int8_t decodeAt(unsigned i, unsigned j) const;

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
	 * Get the alignment string at given pos
	 * @param j  CS position
	 * @return  the alignment at the j-th CS pos
	 */
	string alignAt(unsigned j) const;

	/**
	 * Get the consensus residual at given pos
	 * @param j  CS position
	 * @return the consensus residual as the most frequent one
	 */
	char CSResidualAt(unsigned j) const;

	/**
	 * Get the identity of aligned residuals at given pos
	 * @param j  CS posotion
	 * @return  identity that is the frection of residuals that match to the consensus one
	 */
	double identityAt(unsigned j) const;

	/**
	 * Get the percent of gaps of at given pos
	 * @param j  CS posotion
	 * @return  frection of gaps at this pos
	 */
	double percentGapAt(unsigned j) const;

	/**
	 * Update the count matrices of this object
	 */
	void updateCounts();

	/**
	 * Save this MSA object to a binary file
	 * @param f  the binary output
	 * @return  true if everything is good after saved
	 */
	bool save(std::ostream& f);

	/* static methods */
	/**
	 * Load an MSA file in given format
	 * @param filename  MSA file name
	 * @param format  MSA file format
	 * @return  a newly constructed MSA pointer
	 * @throws invalid_argument if alphabet or format is not known
	 */
	static MSA* load(std::istream& f);

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

	virtual ~MSA() {
		clear();
	}

private:
	/* constructors */
	/**
	 * construct an MSA with given alphabet
	 * @throw invalid_argument if the alphabet is not known
	 */
	explicit MSA(const string& alphabet = "dna") : alphabet(alphabet), abc(SeqCommons::getAlphabetByName(alphabet)),
		numSeq(0), csLen(0), totalNumGap(0), isPruned(false), resCount(NULL), gapCount(NULL) {  }

	/* Disable copy and assignment constructors */
	MSA(const MSA& other);
	MSA& operator=(const MSA& other);

	/* Clear the heap memories */
	void clear();

	/* calculate the CS if not provided by the MSA file */
	void calculateCS();

	string alphabet;
	const DegenAlphabet* abc;
	unsigned numSeq; /* number of sequences */
	unsigned csLen;  /* consensus seq length */
	unsigned long totalNumGap;
	//vector<string> ids;
	string concatMSA; // concatenated MSA
	string CS;        // Consensus Sequence
	bool isPruned; // flag for whether this MS is pruned
	unsigned **resCount; /* Residual count matrix w/ csLen * alphabet-size dimension */
	unsigned *gapCount; /* gap count array w/ CSLen length */
};

inline unsigned long MSA::getMSALen() const {
	return static_cast<unsigned long> (numSeq) * static_cast<unsigned long> (csLen);
}

inline unsigned long MSA::getMSANonGapLen() const {
	return getMSALen() - totalNumGap;
}

inline char MSA::residualAt(unsigned i, unsigned j) const {
/*	if(!(i >= 0 && i < numSeq && j >= 0 && j < csLen))
		throw out_of_range("residual is out of range");*/
	return concatMSA[i * csLen + j];
}

inline int8_t MSA::decodeAt(unsigned i, unsigned j) const {
	return abc->decode(residualAt(i, j));
}

inline string MSA::seqAt(unsigned i) const {
	return concatMSA.substr(i * csLen, csLen);
}

inline string MSA::alignAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // only check j range once
		throw out_of_range("CS pos is out of range");
	string aln;
	for(unsigned i = 0; i < numSeq; ++i)
		aln.push_back(concatMSA[i * csLen + j]);
	return aln;
}

inline MSA* MSA::loadMSAFile(const string& alphabet,
		const string& filename, const string& format) {
	MSA* msa = NULL;
	if(format == "fasta")
		msa = loadFastaFile(alphabet, filename);
	else
		throw invalid_argument("Unsupported MSA file format '" + format + "'");
	/* construct the CS, if it is not set bythe MSA file yet */
	msa->calculateCS();
	return msa;
}

inline bool MSA::saveMSAFile(const string& filename, const string& format) {
	if(format == "fasta")
		return MSA::saveFastaFile(filename);
	else throw invalid_argument("Cannot save MSA to file, unsupported MSA file format " + format);
}

} /* namespace EGriceLab */

#endif /* MSA_H_ */
