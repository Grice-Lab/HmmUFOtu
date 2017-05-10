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
#include "PrimarySeq.h"

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

};

} /* namespace EGriceLab */

#endif /* SRC_SEQUTILS_H_ */
