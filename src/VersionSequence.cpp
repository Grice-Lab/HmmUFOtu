/*
 * ProgVer.cpp
 *
 *  Created on: Oct 4, 2017
 *      Author: zhengqi
 */

#include "VersionSequence.h"
#include "StringUtils.h"

#include <cassert>

namespace EGriceLab {

using namespace std;

ostream& VersionSequence::save(ostream& out) const {
	return StringUtils::saveString(toString(), out, true);
}

istream& VersionSequence::load(istream& in) {
	string verStr;
	StringUtils::loadString(verStr, in);
	parseString(verStr, *this);
	return in;
}

istream& operator>>(istream& in, VersionSequence& ver) {
	string verStr;
	in >> verStr;
	VersionSequence::parseString(verStr, ver);
	return in;
}

void VersionSequence::parseString(const string& str, VersionSequence& ver) {
	sscanf(str.c_str(), "v%d.%d.%d", &ver.majorVer, &ver.minorVer, &ver.buildVer);
}

} /* namespace EGriceLab */
