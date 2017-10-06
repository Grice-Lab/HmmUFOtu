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
 * VersionSequence.cpp
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
