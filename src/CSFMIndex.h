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
#include "divsufsort.h"
#include "WaveletTreeNoptrs.h"

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
	vector<uint16_t> locate(const string& pattern) const;

	/**
	 * Locate the consensus sequence positions of given pattern
	 * @param pattern  the un-coded pattern
	 * @return  a random CS position
	 */
	uint16_t locateOne(const string& pattern) const;

	/**
	 * Locate the consensus sequence positions of given pattern
	 * @param pattern  the un-coded pattern
	 * @return  the first CS position on SA order
	 */
	uint16_t locateFirst(const string& pattern) const;

	static const unsigned RRR_SAMPLE_RATE = 8; /* usa a dense RRR sample rate for maximum efficiency */

private:
	/* private default constructor */
	CSFMIndex();
	/* disable the copy and assignment constructor */
	CSFMIndex(const CSFMIndex& other);
	CSFMIndex& operator=(const CSFMIndex& other);

	const DegenAlphabet* abc;
	uint16_t csLen; /* consensus length */
	//uint8_t* concatSeq; /* concatenated alphabet-encoded non-Gap seq */
	int32_t concatLen; /* total length of concatenated encoded non-gap seq */
	int32_t C[UINT8_MAX + 1];

	string csSeq; /* 1-based consensus seq with dummy position at 0 */
	double* csIdentity; /* 1-based consensus identity */
	uint16_t* csSA; /* 1-based index for mapping position from SA to consensus-seq */
	cds_static::WaveletTreeNoptrs* RRR_bwt; /* Wavelet-Tree transformed BWT string */
};

} /* namespace EGriceLab */

#endif /* CSFMINDEX_H_ */
