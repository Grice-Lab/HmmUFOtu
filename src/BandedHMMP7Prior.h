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
	/* member fields */
	DirichletMixture dmME; /* mixture for match emissions */
	DirichletDensity dmIE; /* density for insertion emissions */
	DirichletDensity dmMT; /* density for match transitions */
	DirichletDensity dmIT; /* density for insertion transitions */
	DirichletDensity dmDT; /* density for deletion transitions */

	/* non-member functions */
	/** read content from input */
	friend istream& operator>>(istream& in, BandedHMMP7Prior& pri);
	/** write content into output */
	friend ostream& operator<<(ostream& out, const BandedHMMP7Prior& pri);
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
	return in;
}


} /* namespace EGriceLab */

#endif /* SRC_BANDEDHMMP7PRIOR_H_ */
