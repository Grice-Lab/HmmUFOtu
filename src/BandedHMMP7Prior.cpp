/*
 * BandedHMMP7Prior.cpp
 *
 *  Created on: Jun 13, 2016
 *      Author: zhengqi
 */

#include "BandedHMMP7Prior.h"

namespace EGriceLab {

ostream& operator<<(ostream& out, const BandedHMMP7Prior& pri) {
	out << "Match emission:" << endl << pri.dmME;
	out << "Insert emission:" << endl << pri.dmIE;
	out << "Match transition:" << endl << pri.dmMT;
	out << "Insert transition:" << endl << pri.dmIT;
	out << "Delete transition:" << endl << pri.dmDT;
	return out;
}

} /* namespace EGriceLab */
