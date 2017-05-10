/*
 * SeqUtils.cpp
 *
 *  Created on: May 10, 2017
 *      Author: zhengqi
 */

#include <cassert>
#include "SeqUtils.h"

namespace EGriceLab {
double SeqUtils::pDist(const DigitalSeq& seq1, const DigitalSeq& seq2,
		DigitalSeq::size_type start, DigitalSeq::size_type end) {
	assert(seq1.getAbc() == seq2.getAbc());

	DigitalSeq::size_type d = 0;
	DigitalSeq::size_type N = 0;
	for(DigitalSeq::size_type i = start; i <= end; ++i) {
		int b1 = seq1[i];
		int b2 = seq2[i];
		if(b1 >= 0 && b2 >= 0) {
			N++;
			if(b1 != b2)
				d++;
		}
	}
	return static_cast<double>(d) / N;
}

} /* namespace EGriceLab */
