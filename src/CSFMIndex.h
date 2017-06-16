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
 * CSFMIndex.h
 *
 *  Created on: Nov 5, 2015
 *      Author: zhengqi
 */

#ifndef CSFMINDEX_H_
#define CSFMINDEX_H_

#include <vector>
#include <fstream>
#include <set>
#include <algorithm>
#include "MSA.h"
#include "BandedHMMCommons.h"
#include "divsufsort.h"
#include "WaveletTreeNoptrs.h"
#include "BitSequence.h"
//#include "Array.h"

namespace EGriceLab {
using std::vector;

/**
 * A Consensus-Sequence FM-index for ultra-fast indexing the consensus positions of a multiple-sequence alignment
 */
class CSFMIndex {
public:
	/* constructors */
	/** Default constructor, zero-initiate all members */
	CSFMIndex() : abc(NULL), gapCh('\0'), csLen(0),
			concatLen(0), C(), csIdentity(NULL), concat2CS(NULL),
			saSampled(NULL), saIdx(NULL), bwt(NULL) {
	}

	/** Virtual destructor */
	virtual ~CSFMIndex() {
		clear();
	}

	/** getters */
	uint16_t getCSLen() const {
		return csLen;
	}

	int32_t getConcatLen() const {
		return concatLen;
	}

	/**
	 * Build an CSFMIndex from a MSA object, old data is removed
	 * @param msa  pointer to an MSA object
	 * @return a fresh allocated CSFMIndex
	 */
	CSFMIndex& build(const MSA& msa);

	/** test whether this CSFMIndex object is fully initiated */
	bool isInitiated() const {
		return abc != NULL && gapCh != '\0' && csLen > 0
				&& concatLen > 0 && C != NULL && concat2CS != NULL
				&& saSampled != NULL && saIdx != NULL && bwt != NULL;
	}

	virtual void clear();

	/**
	 * Save this object to ofstream
	 * @param f  ofstream to save to
	 * @return true if saved successfully
	 */
	std::ofstream& save(std::ofstream& out) const;

	/**
	 * Create and load an CSFMIndex from a stream
	 * @param f  ifstream to load from
	 * @return pointer to fresh created CSFMIndex
	 */
	std::ifstream& load(std::ifstream& in);

	/**
	 * Count how many times a given pattern occurs
	 * @param pattern  the un-coded pattern
	 * @return  number of matches found
	 */
	int32_t count(const string& pattern) const;

	/**
	 * Locate the consensus sequence positions of given pattern
	 * @param pattern  the un-coded pattern
	 * @return  a vector of 1-based consensus locations
	 */
	vector<CSLoc> locate(const string& pattern) const;

	/**
	 * Locate the consensus sequence positions of given pattern
	 * @param pattern  the un-coded pattern
	 * @return  a random CS position
	 */
	CSLoc locateOne(const string& pattern) const;

	/**
	 * Locate the consensus sequence positions of given pattern
	 * @param pattern  the un-coded pattern
	 * @return  the first CS position on SA order
	 */
	CSLoc locateFirst(const string& pattern) const;

	/**
	 * Locate the index of the original sequences (0 .. (concatLen / (csLen + 1)) in the concatSeq of given pattern
	 * @param pattern  the un-coded pattern
	 * @return  a vector of the 0-based indices in which sequences the pattern can be found
	 */
	set<unsigned> locateIndex(const string& pattern) const;

	static const unsigned SA_SAMPLE_RATE = 4;  /* sample rate for SA */
	static const unsigned RRR_SAMPLE_RATE = 8; /* RRR sample rate for BWT */
	static const char sepCh = '\0';

	/* friend functions */
	friend void swap(CSFMIndex& lhs, CSFMIndex& rhs);

private:
	/* disable the copy and assignment constructor */
	CSFMIndex(const CSFMIndex& other);
	CSFMIndex& operator=(const CSFMIndex& other);

	/* private functions */
	/*
	 * Access a given SA loc, either by directly searching the stored value or the next sampled value
	 * @param i  1-index on SA
	 * @return  1-index on concatSeq
	 */
	uint32_t accessSA(uint32_t i) const;

	/**
	 * Extract consensus sequence of given region at concatSeq location
	 * @param start  start on BWT string
	 * @param len  length of BWT string
	 * @return the CS of this region, with gaps filled with default gap characters
	 */
	string extractCS(int32_t start, const string& pattern) const;

	/**
	 * build basic information
	 */
	 void buildBasic(const MSA& msa);

	/**
	 * build a concatSeq and update concat2CS index from a MSA
	 */
	uint8_t* buildConcatSeq(const MSA& msa);

	/** build saSampled, saIdx and BWT from other members */
	void buildBWT(const uint8_t* concatSeq);

	const DegenAlphabet* abc;
	char gapCh;
	uint16_t csLen; /* consensus length */
	//uint8_t* concatSeq; /* concatenated alphabet-encoded non-Gap seq */
	int32_t concatLen; /* total length of concatenated encoded non-gap seq, plus null separators between each individual seq */
	int32_t C[UINT8_MAX + 1]; /* cumulative count of each alphabet frequency, with C[0] as dummy position */

	string csSeq; /* 1-based consensus seq with dummy position at 0 */
	double* csIdentity; /* 1-based consensus identity index */

	uint16_t* concat2CS; /* 1-based concatSeq pos to CS pos, 0 for gap pos on CS */
	uint32_t* saSampled; /* 1-based sampled SA of concatSeq */
	cds_static::BitSequence* saIdx; /* 1-based bit index for telling whether this SA position is sampled */
	cds_static::WaveletTreeNoptrs* bwt; /* Wavelet-Tree transformed BWT string for forward concatSeq */
};

inline void CSFMIndex::clear() {
	//delete[] concatSeq;
	delete[] csIdentity;
	delete[] concat2CS;
	delete[] saSampled;
	delete saIdx;
	delete bwt;
}

} /* namespace EGriceLab */

#endif /* CSFMINDEX_H_ */
