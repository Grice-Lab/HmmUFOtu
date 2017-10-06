/*
 * ProgVer.cpp
 *
 *  Created on: Oct 4, 2017
 *      Author: zhengqi
 */

#include "VersionSequence.h"

#include <cassert>

namespace EGriceLab {

using namespace std;

ostream& VersionSequence::save(ostream& out) const {
	out.write((const char*) &majorVer, sizeof(int));
	out.write((const char*) &minorVer, sizeof(int));
	out.write((const char*) &buildVer, sizeof(int));
	return out;
}

istream& VersionSequence::load(istream& in) {
	in.read((char*) &majorVer, sizeof(int));
	in.read((char*) &minorVer, sizeof(int));
	in.read((char*) &buildVer, sizeof(int));
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
