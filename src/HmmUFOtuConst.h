/*
 * HmmUFOtuConst.h
 *
 *  Created on: Jun 3, 2015
 *      Author: zhengqi
 */
#include<string>
#include<iostream>
#include<limits>
#include<cassert>

#ifndef HMMUFOTUCONST_H_
#define HMMUFOTUCONST_H_

namespace EGriceLab {
using std::string;
using std::istream;
using std::ostream;

/* includes for some program-wide constants and definitions */
const std::string progName = "HmmUFOtu";
const std::string progVersion = "v1.01";

const double inf = std::numeric_limits<double>::infinity();
const double infV = -inf;

/* stand-alone functions */
/**
 * Compare two versions
 */
int cmpVersion(const string& ver1, const string& ver2);

/**
 * Read progName from an input stream
 */
istream& readProgName(istream& in, string& name);

/**
 * Read progVersion from an input stream
 */
istream& readProgVersion(istream& in, string& version);

/**
 * Write progName to an output stream
 */
ostream& writeProgName(ostream& out, const string& name);

/**
 * Write progVersion to an output stream
 */
ostream& writeProgVersion(ostream& out, const string& version);

} /* end namespace EGriceLab */
#endif /* HMMUFOTUCONST_H_ */
