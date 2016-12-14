/*
 * BandedHMMP7Bg.h
 *
 *  Created on: May 11, 2015
 *      Author: zhengqi
 */

#ifndef BANDEDHMMP7BG_H_
#define BANDEDHMMP7BG_H_
#include <cmath>
#include <cassert>
#include <limits>
#include <eigen3/Eigen/Dense>
#include "IUPACNucl.h"
#include "SeqCommons.h"

namespace EGriceLab {
using Eigen::Vector4d;
/*
 * A class to represent the background transition and emission distributions of a null banded p7 model,
 * consisting the bg state G and dummy end state F, like the T in plan 7 architecture
 */
class BandedHMMP7Bg {
public:
	/* constructors */
	/*
	 * constructor with given size
	 */
	explicit BandedHMMP7Bg(int size, const IUPACNucl* abc = SeqCommons::nuclAbc)
	: K(size), nuclAbc(abc) {
		init_bgFreq();
		init_transPr();
	}
	/* member methods */
	/**
	 * return the background transition prob between G states
	 */
	float getBgTransPr() const {
		return p1;
	}

	/**
	 * return the background transition lods between G states
	 */
	float getBgTransLods() const {
		return log(p1);
	}

	/*
	 * return the background emission vector at G state
	 */
	Vector4d getBgEmitPr() const {
		return bgFreq;
	}

	/*
	 * return the background emission vector at G state in log scale
	 */
	Vector4d getBgEmitLogPr() const {
		return bgFreq.array().log();
	}

	/*
	 * return the background emission probability of Aphabet i
	 */
	float getBgEmitPr(int i) const {
		return bgFreq(i);
	}

	/**
	 * reset the size of this background model, adjusting the transition prob accordingly
	 */
	void setSize(int size);
	/**
	 * set the background nucleotide frequencies using observed frequencies or count
	 * @param freq  the observed frequencies or count of each nucleotide
	 */
	void setBgFreq(const Vector4d& q);

private:
	/* private member functions */
	void init_bgFreq();
	void init_transPr();
	//void init_emisPr();

	int K; // profile size
	const DegenAlphabet* nuclAbc;
	Vector4d bgFreq; // null background frequencies of each nuclotide bases
	double p1; // null transition distribution of G->G, which is 1 - p0 = 1 - transBG

	static const int MIN_BG_K = 350; /* min profile length used to set bg transition probability */
};

} /* namespace EGriceLab */

#endif /* BANDEDHMMP7BG_H_ */
