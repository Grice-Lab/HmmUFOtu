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
#include "MSA.h"
#include "BandedHMMCommons.h"
#include "divsufsort.h"
#include "WaveletTreeNoptrs.h"
#include "Array.h"

namespace EGriceLab {

/**
 * A Consensus-Sequence FM-index for ultra-fast indexing the consensus positions of a multiple-sequence alignment
 */
class CSFMIndex {
public:
	virtual ~CSFMIndex();

	/**
	 * Build an CSFMIndex from a MSA object
	 * @param msa  pointer to an MSA object
	 * @return a fresh allocated CSFMIndex
	 */
	static CSFMIndex* build(const MSA* msa);

	/**
	 * Save this object to ofstream
	 * @param f  ofstream to save to
	 * @return true if saved successfully
	 */
	bool save(std::ofstream& f) const;

	/**
	 * Create and load an CSFMIndex from a stream
	 * @param f  ifstream to load from
	 * @return pointer to fresh created CSFMIndex
	 */
	static CSFMIndex* load(std::ifstream& f);

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

	static const unsigned RRR_SAMPLE_RATE = 8; /* usa a dense RRR sample rate for maximum efficiency */

private:
	/* private default constructor */
	CSFMIndex();
	/* disable the copy and assignment constructor */
	CSFMIndex(const CSFMIndex& other);
	CSFMIndex& operator=(const CSFMIndex& other);

	/* private functions */
	/**
	 * Extract consensus sequence of given region at concatSeq location
	 * @param start  start on BWT string
	 * @param len  length of BWT string
	 * @return the CS of this region, with gaps filled with default gap characters
	 */
	string extractCS(int32_t start, int32_t len) const;

	/**
	 * Transfer loc to reverse orientation on concatSeq
	 * @param loc location on this current direction
	 * return location on the reverse direction
	 */
	int32_t reverseLoc(int32_t loc);

	const DegenAlphabet* abc;
	char gapCh;
	uint16_t csLen; /* consensus length */
	//uint8_t* concatSeq; /* concatenated alphabet-encoded non-Gap seq */
	int32_t concatLen; /* total length of concatenated encoded non-gap seq */
	int32_t C[UINT8_MAX + 1]; /* cumulative count of each alphabet frequency, with C[0] as dummy position */

	string csSeq; /* 1-based consensus seq with dummy position at 0 */
	double* csIdentity; /* 1-based consensus identity */
	uint16_t* csSA; /* 1-based index for mapping position from SA to consensus-seq */
	cds_static::WaveletTreeNoptrs* fwt_bwt; /* Wavelet-Tree transformed BWT string for forward concatSeq */
	cds_static::WaveletTreeNoptrs* rev_bwt; /* Wavelet-Tree transformed BWT string for reverse concatSeq */
};

inline int32_t CSFMIndex::reverseLoc(int32_t loc) {
	return concatLen - 1 - loc;
}

} /* namespace EGriceLab */

#endif /* CSFMINDEX_H_ */
