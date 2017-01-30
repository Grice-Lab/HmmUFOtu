/*
 * BandedHMMP7Prior.cpp
 *
 *  Created on: Jun 13, 2016
 *      Author: zhengqi
 */

#include "BandedHMMP7Prior.h"

namespace EGriceLab {

const double BandedHMMP7Prior::DEFAULT_ABS_EPS_COST = 0;
const double BandedHMMP7Prior::DEFAULT_REL_EPS_COST = 1e-6;
const double BandedHMMP7Prior::DEFAULT_ABS_EPS_PARAMS = 0;
const double BandedHMMP7Prior::DEFAULT_REL_EPS_PARAMS = 1e-6;

ostream& operator<<(ostream& out, const BandedHMMP7Prior& pri) {
	out << "Match emission:" << endl << pri.dmME;
	out << "Insert emission:" << endl << pri.dmIE;
	out << "Match transition:" << endl << pri.dmMT;
	out << "Insert transition:" << endl << pri.dmIT;
	out << "Delete transition:" << endl << pri.dmDT;
	return out;
}

} /* namespace EGriceLab */
