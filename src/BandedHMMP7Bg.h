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
#include <Eigen/Dense>

#include "AlphabetFactory.h"
#include "IUPACNucl.h"

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
	explicit BandedHMMP7Bg(int size, const DegenAlphabet* abc = AlphabetFactory::getAlphabetByName("DNA"))
	: K(size), nuclAbc(abc) {
		init_bgFreq();
		init_transPr();
	}
	/* member methods */
	/**
	 * return the background transition prob between G states
	 */
	double getBgTransPr() const {
		return p1;
	}

	/** return the background termination prob between G and F state */
	double getBgTermPr() const {
		return 1 - p1;
	}

	/**
	 * return the background transition lods between G states
	 */
	double getBgTransLods() const {
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
	double getBgEmitPr(int i) const {
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
