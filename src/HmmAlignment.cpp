/*
 * HmmAlignment.cpp
 *
 *  Created on: Feb 1, 2018
 *      Author: zhengqi
 */

#include "HmmAlignment.h"
#include "BandedHMMP7.h"

namespace EGriceLab {
namespace HmmUFOtu {

const string HmmAlignment::TSV_HEADER = "K\tL\tseq_start\tseq_end\thmm_start\thmm_end\tCS_start\tCS_end\tcost\talignment";

HmmAlignment& HmmAlignment::merge(const HmmAlignment& other) {
	if(isCompatitable(other)) {
		/* merge seq loc */
		if(other.seqStart < seqStart)
			seqStart = other.seqStart;
		if(other.seqEnd > seqEnd)
			seqEnd = other.seqEnd;
		/* merge HMM loc */
		if(other.hmmStart < hmmStart)
			hmmStart = other.hmmStart;
		if(other.hmmEnd > hmmEnd)
			hmmEnd = other.hmmEnd;
		/* merge CS loc */
		if(other.csStart < csStart)
			csStart = other.csStart;
		if(other.csEnd > csEnd)
			csEnd = other.csEnd;
		/* add cost */
		cost += other.cost;
		/* merge aligned seq */
		for(string::size_type i = 0; i < L; ++i)
			if(align[i] == BandedHMMP7::PAD_SYM && other.align[i] != BandedHMMP7::PAD_SYM) /* this align has priority */
				align[i] = other.align[i];
	}
	return *this;
}

ostream& operator<<(ostream& out, const HmmAlignment& hmmAln) {
	out << hmmAln.K << "\t" << hmmAln.L << "\t" <<
			hmmAln.seqStart << "\t" << hmmAln.seqEnd << "\t" <<
			hmmAln.hmmStart << "\t" << hmmAln.hmmEnd << "\t" <<
			hmmAln.csStart << "\t" << hmmAln.csEnd << "\t" <<
			hmmAln.cost << "\t" << hmmAln.align;
	return out;
}

istream& operator>>(istream& in, HmmAlignment& hmmAln) {
	in >> hmmAln.K >> hmmAln.L >>
	hmmAln.seqStart >> hmmAln.seqEnd >>
	hmmAln.hmmStart >> hmmAln.hmmEnd >>
	hmmAln.csStart >> hmmAln.csEnd >>
	hmmAln.cost >> hmmAln.align;
	return in;
}


} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
