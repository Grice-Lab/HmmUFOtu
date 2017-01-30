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
#include "StringUtils.h"
#include "DirichletModel.h"
#include "DirichletDensity.h"
#include "DirichletMixture.h"

namespace EGriceLab {
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
		setAbsEpsParams(DEFAULT_ABS_EPS_COST);
		setRelEpsParams(DEFAULT_REL_EPS_COST);
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

inline istream& operator>>(istream& in, BandedHMMP7Prior& pri) {
	string head;
	while(std::getline(in, head)) {
		if(StringUtils::startsWith(head, "Match emission:"))
			in >> pri.dmME;
		else if(StringUtils::startsWith(head, "Insert emission:"))
			in >> pri.dmIE;
		else if(StringUtils::startsWith(head, "Match transition:"))
			in >> pri.dmMT;
		else if(StringUtils::startsWith(head, "Insert transition:"))
			in >> pri.dmIT;
		else if(StringUtils::startsWith(head, "Delete transition:"))
			in >> pri.dmDT;
		else
			continue;
	}
	if(!(pri.dmME.getK() > 0 && pri.dmIE.getK() > 0
			&& pri.dmMT.getK() > 0 && pri.dmIT.getK() > 0 && pri.dmDT.getK() > 0)) {
		std::cerr << "Empty or partial BandedHMMP7Prior input" << endl;
		in.setstate(std::ios_base::failbit);
	}
	return in;
}

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

} /* namespace EGriceLab */

#endif /* SRC_BANDEDHMMP7PRIOR_H_ */
