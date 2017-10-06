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
 * HmmUFOtuEnv.cpp
 *
 *  Created on: Jul 15, 2016
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include <cstring>
#include "HmmUFOtuEnv.h"
#include "HmmUFOtuConst.h"
#include "StringUtils.h"

namespace EGriceLab {
using namespace std;

int VERBOSE_LEVEL = LOG_WARNING; /* DEFAULT VERBOSE LEVEL */
const VersionSequence progVer("v1.2.2");
const string projectURL = "https://github.com/Grice-Lab/HmmUFOtu";

void printVersion(const string& app, ostream& out) {
	out << app << ": " << progVer << std::endl;
	out << "Package: " << progName << " " << progVer << std::endl;
}

ostream& saveProgInfo(ostream& out) {
	StringUtils::saveString(progName, out, progName.length()); /* save name with known length */
	progVer.save(out); /* save version */
	return out;
}

istream& loadProgInfo(istream& in) {
	/* load program info */
	string pname;
	VersionSequence pver;
	StringUtils::loadString(pname, in, progName.length());

	/* load name */
	if(in.bad()) {
		cerr << "Unable to load database file: " << ::strerror(errno) << endl;
		return in;
	}
	if(progName != pname) {
		cerr << "Not an valid database file of " << progName << endl;
		in.setstate(ios_base::badbit);
		return in;
	}

	/* load version */
	pver.load(in);
	if(in.bad()) {
		cerr << "Unrecognized " << progName << " version: " << ::strerror(errno) << endl;
		return in;
	}

	if(!(progVer >= pver)) {
		cerr << "You are using an old version of " << getProgFullName(progName, progVer)
				<< " to read a newer database file that is build by " << getProgFullName(pname, pver)
				<< " please download the latest program from '"
				<< projectURL << "'" << endl;
		in.setstate(ios_base::badbit);
		return in;
	}

	return in;
}

ostream& writeProgInfo(ostream& out, const string& info) {
	out << "# " << progName << " " << progVer << info << endl;
	return out;
}

istream& readProgInfo(istream& in) {
	string header;
	std::getline(in, header);
	/* check program info */
	char pname[MAX_NAME_LENGTH], ver[MAX_NAME_LENGTH];
	if(sscanf(header.c_str(), "# %s %s", pname, ver) != 2) {
		cerr << "Unrecognized input file for " << progName << endl;
		in.setstate(ios_base::badbit);
		return in;
	}

	if(progName != pname) {
		cerr << "Not an valid input file of " << progName << endl;
		in.setstate(ios_base::badbit);
		return in;
	}

	VersionSequence pver(ver);
	if(!(progVer >= pver)) {
		cerr << "You are using an old version of " << getProgFullName(progName, progVer)
				<< " to read a newer input file that is build by " << getProgFullName(pname, pver)
				<< " please download the latest program from '"
				<< projectURL << "'" << endl;
		in.setstate(ios_base::badbit);
		return in;
	}

	return in;
}

}
