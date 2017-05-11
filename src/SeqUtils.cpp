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
	assert(seq1.length() == seq2.length());

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

double SeqUtils::pDist(const string& seq1, const string& seq2,
		string::size_type start, string::size_type end) {
	assert(seq1.length() == seq2.length());
	string::size_type d = 0;
	for(string::size_type i = start; i <= end; ++i)
		if(seq1[i] != seq2[i])
			d++;
	return static_cast<double>(d) / (end - start + 1);
}

double SeqUtils::pDist(const string& seq1, const string& seq2,
		const DegenAlphabet* abc, string::size_type start,
		string::size_type end) {
	assert(seq1.length() == seq2.length());
	string::size_type d = 0;
	for(string::size_type i = start; i <= end; ++i)
		if(!abc->isMatch(seq1[i], seq2[i]))
			d++;
	return static_cast<double>(d) / (end - start + 1);
}

double SeqUtils::pDist(const string& seq1, const DigitalSeq& seq2, size_t start,
		size_t end) {
	assert(seq1.length() == seq2.length());
	const DegenAlphabet* abc = seq2.getAbc();
	size_t d = 0;
	for(size_t i = start; i <= end; ++i)
		if(!abc->isMatch(seq1[i], seq2[i]))
			d++;
	return static_cast<double>(d) / (end - start + 1);
}

} /* namespace EGriceLab */
