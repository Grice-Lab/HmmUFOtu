/*
 * HmmAlignment.h
 *
 *  Created on: Feb 1, 2018
 *      Author: zhengqi
 */

#ifndef SRC_HMMALIGNMENT_H_
#define SRC_HMMALIGNMENT_H_
#include <string>
#include <iostream>
#include "HmmUFOtuConst.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::string;
using std::istream;
using std::ostream;

struct HmmAlignment {
	/* constructors */
	/** default constructor */
	HmmAlignment() {  }

	/** construct from given data */
	HmmAlignment(int K, int L,
			int seqStart, int seqEnd, int hmmStart, int hmmEnd,
			int csStart, int csEnd, double cost, const string& align)
	: K(K), L(L),
	  seqStart(seqStart), seqEnd(seqEnd), hmmStart(hmmStart), hmmEnd(hmmEnd),
	  csStart(csStart), csEnd(csEnd), cost(cost), align(align)
	{  }

	/* member methods */
	bool isValid() const {
		return 0 < seqStart && seqStart <= seqEnd &&
				0 < hmmStart && hmmStart <= hmmEnd && hmmEnd <= K &&
				0 < csStart && csStart <= csEnd && csEnd <= L &&
				cost >= 0 && cost != inf && L == align.length();
	}

	bool isCompatitable(const HmmAlignment& other) const {
		return K == other.K && L == other.L;
	}

	/**
	 * Merge this alignment with another alignment, or do nothing if they are not compatitable
	 */
	HmmAlignment& merge(const HmmAlignment& otherAln);

	/* static methods */
	/**
	 * Merge two HmmAlignments
	 * @return the merged alignment if compatitable, or a copy of the first alignment if not
	 */
	static HmmAlignment merge(const HmmAlignment& aln1, const HmmAlignment& aln2);

	/* non-member friend functions */
	/** write to a text output */
	friend ostream& operator<<(ostream& out, const HmmAlignment& hmmAln);

	/** read from a text input */
	friend istream& operator>>(istream& in, HmmAlignment& hmmAln);

	/* member fields */
	int K; /* HMM profile size */
	int L; /* concensus size */
	int seqStart, seqEnd; /* 1-based seq coordinates */
	int hmmStart, hmmEnd; /* 1-based HMM profile coordinates */
	int csStart, csEnd; /* 1-based consensus coordinates */
	double cost; /* HMM align cost */
	string align; /* alignmented seq */

	/* static fields */
	static const string TSV_HEADER;
};

inline HmmAlignment HmmAlignment::merge(const HmmAlignment& aln1, const HmmAlignment& aln2) {
	HmmAlignment alnMerged(aln1); /* make a local copy */
	return alnMerged.merge(aln2);
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_HMMALIGNMENT_H_ */
