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
#include <Eigen/Dense>
#include <cmath>

#include "AlphabetFactory.h"
#include "StringUtils.h"
#include "PrimarySeq.h"
#include "DigitalSeq.h"
#include "ProgLog.h"

namespace EGriceLab {
namespace HmmUFOtu {

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

	const string& getSeqName(unsigned i) const {
		return seqNames[i];
	}

	bool pruned() const {
		return isPruned;
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
	 * @throws out_of_range exception if out of range
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
	 * @throws  out_of_range exception of j is out of range
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
	 * Get the weighted fraction of gaps of at given pos
	 * @param j  CS position
	 * @return  weighted fraction of gaps at this pos
	 */
	double gapWFrac(unsigned j) const;

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
	 * Get the residual frequency of each base of the entire MSA
	 */
	VectorXd resFreq() const;

	/**
	 * Get the weighted residual frequency of each base of the entire MSA
	 */
	VectorXd resWFreq() const;

	/**
	 * Update the count matrices of this object, also update auxiliary indices
	 */
	void updateRawCounts();

	/**
	 * Update the seqWeight of this object using the Henikoffs' algorithm (1994)
	 */
	void updateSeqWeight();

	/**
	 * Update the count matrices of this object
	 */
	void updateWeightedCounts();

	/**
	 * Scale every seq weight by a constant
	 */
	void sclaleWeight(double r);

	/**
	 * Save raw object data to output
	 */
	std::ostream& save(std::ostream& out) const;

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

	/**
	 * Load raw object data from input
	 */
	std::istream& load(std::istream& in);

	/**
	 * Load an MSA binary file
	 * @param abc  alphabet of input
	 * @param in  input stream
	 * @param format  input format
	 * @return  number of MSA sequences successfully read in
	 */
	long loadMSA(const DegenAlphabet* abc, istream& in, const string& format);

	/**
	 * Load an MSA binary file
	 * @param alphabet  alphabet name of input
	 * @param in  input stream
	 * @param format  input format
	 * @return  number of MSA sequences successfully read in
	 */
	long loadMSA(const string& alphabet, istream& in, const string& format) {
		return loadMSA(AlphabetFactory::getAlphabetByName(alphabet), in, format);
	}

	/**
	 * Load an MSA file in fasta format
	 * @param filename  MSA file name
	 * @return  a newly constructed MSA pointer
	 * @throws invalid_argument if alphabet or format is not known
	 */
	long loadMSAFasta(const DegenAlphabet* abc, istream& in);

	/* constructors */
	/**
	 * construct an MSA with given alphabet
	 * @throw invalid_argument if the alphabet is not known
	 */
	explicit MSA(const string& alphabet = "dna") : alphabet(alphabet), abc(AlphabetFactory::getAlphabetByName(alphabet)),
		numSeq(0), csLen(0), isPruned(false)
	{  }

	/* Clear the heap memories */
	void clear();

	/* reset count values */
	void resetRawCount();

	/* reset seq weights */
	void resetSeqWeight();

	/* reset weighted count values */
	void resetWeightedCount();

	/* calculate the CS if not provided by the MSA file */
	void calculateCS();

private:
	string alphabet;
	const DegenAlphabet* abc; /* stored abc const pointer that guarenteed to be a global variable */
	string name;
	unsigned numSeq; /* number of sequences */
	unsigned csLen;  /* consensus seq length */
	vector<string> seqNames; /* seq names stored in their occurring order */
	//vector<string> ids;
	string concatMSA; // concatenated MSA
	string CS;        // Consensus Sequence
	bool isPruned; // flag for whether this MS is pruned
	/* auxiliary data to remember each sequence start, end and length (non-gapped) */
	vector<int> startIdx; /* start position on CS */
	vector<int> endIdx; /* end position on CS */
	vector<int> lenIdx; /* unmapped length */

	/* matrix/vector for residual & gap count */
	MatrixXi resCount; /* Residual count matrix w/ alphabet-size X CSLen dimension */
	VectorXi gapCount; /* gap count array w/ CSLen length */
	VectorXd seqWeight; /* Sequence weight for each seq w/ numSeq length */
	MatrixXd resWCount; /* weighted residual count matrix w/ alphabet-size X CSLen dimension */
	VectorXd gapWCount; /* weighted gap count array w/ CSLen length */

	/* static members */
public:
	static const double DEFAULT_CONSENSUS_FRAC;
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
	return concatMSA.at(i * csLen + j);
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
	string aln;
	for(unsigned i = 0; i < numSeq; ++i)
		aln.push_back(concatMSA.at(i * csLen + j));
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

inline VectorXd MSA::resFreq() const {
	VectorXd freq = resCount.rowwise().sum().cast<double>();
	return freq / freq.sum();
}

inline VectorXd MSA::resWFreq() const {
	VectorXd freq = resWCount.rowwise().sum();
	return freq / freq.sum();
}

inline void MSA::sclaleWeight(double r) {
	seqWeight *= r;
	updateWeightedCounts();
}

inline long MSA::loadMSA(const DegenAlphabet* abc,
		istream& in, const string& format) {
	if(format == "fasta")
		return loadMSAFasta(abc, in);
	else {
		errorLog << "Unsupported MSA file format '" + format + "'";
		return -1;
	}
}

inline bool MSA::saveMSAFile(const string& filename, const string& format) {
	if(format == "fasta")
		return MSA::saveFastaFile(filename);
	else throw invalid_argument("Cannot save MSA to file, unsupported MSA file format " + format);
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* MSA_H_ */
