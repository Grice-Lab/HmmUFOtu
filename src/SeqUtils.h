/*
 * SeqUtils.h
 *  utility functions provided for Sequences
 *  Created on: May 10, 2017
 *      Author: zhengqi
 */

#ifndef SRC_SEQUTILS_H_
#define SRC_SEQUTILS_H_
#include <string>
#include "DigitalSeq.h"

namespace EGriceLab {

class SeqUtils {
public:
	/* static methods */
	/**
	 * calculate the p-distance between two aligned DigitalSeq in given region [start, end]
	 */
	static double pDist(const DigitalSeq& seq1, const DigitalSeq& seq2,
			DigitalSeq::size_type start, DigitalSeq::size_type end);

	/** calculate the p-distance between two aligned DigitalSeq */
	static double pDist(const DigitalSeq& seq1, const DigitalSeq& seq2) {
		return pDist(seq1, seq2, 0, seq1.length() - 1);
	}

	/** calculate the p-distance between two strings in a given region [start, end] */
	static double pDist(const string& seq1, const string& seq2,
			string::size_type start, string::size_type end);

	/** calculate the p-distance between two strings */
	static double pDist(const string& seq1, const string& seq2) {
		return pDist(seq1, seq2, 0, seq1.length() - 1);
	}

	/** calculate the p-distance between two strings in a given region [start, end], allow Degenerated characters */
	static double pDist(const string& seq1, const string& seq2, const DegenAlphabet* abc,
			string::size_type start, string::size_type end);

	/** calculate the p-distance between two strings */
	static double pDist(const string& seq1, const string& seq2, const DegenAlphabet* abc) {
		return pDist(seq1, seq2, abc, 0, seq1.length() - 1);
	}

	/** calculate the p-distance between two strings in a given region [start, end], allow Degenerated characters */
	static double pDist(const string& seq1, const DigitalSeq& seq2,
			size_t start, size_t end);

	/** calculate the p-distance between two strings */
	static double pDist(const string& seq1, const DigitalSeq& seq2) {
		return pDist(seq1, seq2, 0, seq1.length() - 1);
	}

};

} /* namespace EGriceLab */

#endif /* SRC_SEQUTILS_H_ */
