/*
 * ProgEnv.h
 *
 *  Created on: May 23, 2018
 *      Author: zhengqi
 */

#ifndef PROGENV_H_
#define PROGENV_H_

#include <string>
#include <iostream>

#include "VersionSequence.h"
#include "ProgLog.h"

namespace EGriceLab {
/* per-application variables */
extern const std::string progName;
extern const VersionSequence progVer;
extern const string projectURL;
const int MAX_PROG_LENGTH = 4096;

/**
 * show program and package version
 */
void printVersion(const string& app, ostream& out = std::cerr);

/**
 * get full program name
 */
inline string getProgFullName(const string& name, const VersionSequence& ver) {
	return name + "-" + ver.toString();
}

/** save program info to a binary output */
ostream& saveProgInfo(ostream& out);

/** load program info from a binary input and check agains known values */
istream& loadProgInfo(istream& in);

/** load program info from a binary input and check agains known values */
istream& loadProgInfo(istream& in, VersionSequence& pver);

/** write progInfo and additional information to a text output */
ostream& writeProgInfo(ostream& out, const string& info = "");

/** read progInfo and additional information from a text input */
istream& readProgInfo(istream& in);

/** read progInfo and additional information from a text input */
istream& readProgInfo(istream& in, VersionSequence& pver);

} /* namespace EGriceLab */

#endif /* PROGENV_H_ */
