/*
 * ProgEnv.cpp
 *
 *  Created on: May 23, 2018
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include <cstring>
#include <cerrno>
#include "ProgEnv.h"
#include "StringUtils.h"

namespace EGriceLab {
using namespace std;
int VERBOSE_LEVEL = LOG_WARNING; /* DEFAULT VERBOSE LEVEL */

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

istream& loadProgInfo(istream& in, VersionSequence& pver) {
	/* load program info */
	string pname;
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
	char pname[MAX_PROG_LENGTH], ver[MAX_PROG_LENGTH];
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

istream& readProgInfo(istream& in, VersionSequence& pver) {
	string header;
	std::getline(in, header);
	/* check program info */
	char pname[MAX_PROG_LENGTH], ver[MAX_PROG_LENGTH];
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

	pver = VersionSequence(ver);
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

} /* namespace EGriceLab */

