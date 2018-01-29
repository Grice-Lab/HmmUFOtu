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
 * BandedHMMP7Prior.h
 *
 *  A POD class storing the Dirichlet priors of a BandedHMMP7 model
 *  Created on: Jun 13, 2016
 *      Author: zhengqi
 */

#ifndef SRC_BANDEDHMMP7PRIOR_H_
#define SRC_BANDEDHMMP7PRIOR_H_

#include <string>
#include <iostream>
#include "HmmUFOtuConst.h"
#include "StringUtils.h"
#include "DirichletModel.h"
#include "DirichletDensity.h"
#include "DirichletMixture.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::istream;
using std::ostream;
using std::endl;
using std::string;
using Math::DirichletModel;
using Math::DirichletDensity;
using Math::DirichletMixture;

struct BandedHMMP7Prior {
	/* constructors */
	BandedHMMP7Prior() {
		setMaxIter(DEFAULT_MAX_ITER);
		setAbsEpsCost(DEFAULT_ABS_EPS_COST);
		setRelEpsCost(DEFAULT_REL_EPS_COST);
		setAbsEpsParams(DEFAULT_ABS_EPS_PARAMS);
		setRelEpsParams(DEFAULT_REL_EPS_PARAMS);
	}

	/* member fields */
	DirichletMixture dmME; /* mixture for match emissions */
	DirichletDensity dmIE; /* density for insertion emissions */
	DirichletDensity dmMT; /* density for match transitions */
	DirichletDensity dmIT; /* density for insertion transitions */
	DirichletDensity dmDT; /* density for deletion transitions */

	/* member functions */
	/* convenient setters to forward calls to underlying models */
	void setDims(int K, int L);

	void setMaxIter(int maxIter);
	void setAbsEpsCost(double eps);
	void setRelEpsCost(double eps);
	void setAbsEpsParams(double eps);
	void setRelEpsParams(double eps);

	/* non-member functions */
	/** read content from input */
	friend istream& operator>>(istream& in, BandedHMMP7Prior& pri);
	/** write content into output */
	friend ostream& operator<<(ostream& out, const BandedHMMP7Prior& pri);

	/* static members */
	static const int DEFAULT_MAX_ITER = 0;
	static const double DEFAULT_ABS_EPS_COST;
	static const double DEFAULT_REL_EPS_COST;
	static const double DEFAULT_ABS_EPS_PARAMS;
	static const double DEFAULT_REL_EPS_PARAMS;
};

inline void BandedHMMP7Prior::setDims(int K, int L) {
	/* set the # of parameters */
	dmME.setDims(K, L);
	dmIE.setK(K);
	dmMT.setK(3); /* M->M, M->I, M-D */
	dmIT.setK(2); /* I->M, I->I */
	dmDT.setK(2); /* D->M, D->D */
}

inline void BandedHMMP7Prior::setMaxIter(int maxIter) {
	dmME.setMaxIter(maxIter);
	dmIE.setMaxIter(maxIter);
	dmMT.setMaxIter(maxIter);
	dmIT.setMaxIter(maxIter);
	dmDT.setMaxIter(maxIter);
}

inline void BandedHMMP7Prior::setAbsEpsCost(double eps) {
	dmME.setAbsEpsCost(eps);
	dmIE.setAbsEpsCost(eps);
	dmMT.setAbsEpsCost(eps);
	dmIT.setAbsEpsCost(eps);
	dmDT.setAbsEpsCost(eps);
}

inline void BandedHMMP7Prior::setRelEpsCost(double eps) {
	dmME.setRelEpsCost(eps);
	dmIE.setRelEpsCost(eps);
	dmMT.setRelEpsCost(eps);
	dmIT.setRelEpsCost(eps);
	dmDT.setRelEpsCost(eps);
}

inline void BandedHMMP7Prior::setAbsEpsParams(double eps) {
	dmME.setAbsEpsParams(eps);
	dmIE.setAbsEpsParams(eps);
	dmMT.setAbsEpsParams(eps);
	dmIT.setAbsEpsParams(eps);
	dmDT.setAbsEpsParams(eps);
}

inline void BandedHMMP7Prior::setRelEpsParams(double eps) {
	dmME.setRelEpsParams(eps);
	dmIE.setRelEpsParams(eps);
	dmMT.setRelEpsParams(eps);
	dmIT.setRelEpsParams(eps);
	dmDT.setRelEpsParams(eps);
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_BANDEDHMMP7PRIOR_H_ */
