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
 * DigitalSeq.cpp
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include "DigitalSeq.h"
#include "StringUtils.h"
#include "AlphabetFactory.h"

using namespace std;

namespace EGriceLab {

DigitalSeq::DigitalSeq(const DegenAlphabet* abc, const string& name, const string& str) :
				abc(abc), name(name) {
	for(string::const_iterator it = str.begin(); it != str.end(); ++it) {
		char c = ::toupper(*it);
		if(abc->isValid(c))
			push_back(abc->encode(c)); // use encoded values
	}
}

DigitalSeq::DigitalSeq(const PrimarySeq& seq) :
	abc(seq.getAbc()), name(seq.getId()) {
	append(seq.getSeq());
}

string DigitalSeq::toString() const {
	string str;
	for(DigitalSeq::const_iterator it = begin(); it != end(); ++it)
		str.push_back(abc->decode(*it));
	return str;
}

DigitalSeq DigitalSeq::revcom() const {
	if(!abc->hasComplement())
		throw std::invalid_argument("Sequence alphabet " + abc->getName() + " does not support reverse-complement");
	DigitalSeq revcomSeq(abc, name); // make an empty copy with same DegebAlphabet and name
	for(DigitalSeq::const_reverse_iterator rit = rbegin(); rit != rend(); ++rit)
		revcomSeq.push_back(abc->encode(abc->getComplementSymbol(abc->decode(*rit))));
	return revcomSeq;
}

string DigitalSeq::join(const string& sep) {
	ostringstream ostr;
	for(const_iterator it = begin(); it != end(); ++it) {
		if(it != begin())
			ostr << sep;
		ostr << *it;
	}
	return ostr.str();
}

DigitalSeq& EGriceLab::DigitalSeq::append(const string& str) {
	for(string::const_iterator it = str.begin(); it != str.end(); ++it) {
		char c = ::toupper(*it);
		if(abc->isValid(c))
			push_back(abc->encode(c));
	}
	return *this;
}

/*
 * compare two DigitalSeq
 * return true if and only if all residuals are equal and are the same Alphabet
 */
bool operator==(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	if(!(lhs.abc == rhs.abc /* same object */ || *lhs.abc == *rhs.abc /* equal object */))
		return false;
	if(lhs.length() != rhs.length())
		return false;
	for(DigitalSeq::size_type i = 0; i != lhs.length(); ++i)
		if(lhs[i] != rhs[i])
			return false;
	return true;
}

/*
 * compare two DigitalSeq strict weak order, based on lexical order of the decoded string
 * return true if and only if lhs is strictly less than rhs
 */
bool operator<(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	return lhs.toString() < rhs.toString();
}


ostream& operator<<(ostream& os, const DigitalSeq& dSeq) {
	for(DigitalSeq::const_iterator it = dSeq.begin(); it != dSeq.end(); ++it)
		os << dSeq.abc->decode(*it);
	return os;
}

ostream& DigitalSeq::save(ostream& out) const {
	/* empty flag */
	bool initiated = abc != NULL;
	out.write((const char*) &initiated, sizeof(bool));
	if(!initiated)
		return out;

	/* save sizes */
	string alphabet = abc->getName();
	string::size_type nAlphabet = alphabet.length();
	string::size_type nName = name.length();
	DigitalSeq::size_type len = length();

	out.write((const char*) &nAlphabet, sizeof(string::size_type));
	out.write((const char*) &nName, sizeof(string::size_type));
	out.write((const char*) &len, sizeof(DigitalSeq::size_type));

	/* save basic info */
	StringUtils::saveString(alphabet, out);
	StringUtils::saveString(name, out);

	/* save seq */
	StringUtils::saveString(*this, out);
	return out;
}

istream& DigitalSeq::load(istream& in) {
	/* read flag */
	bool initiated;
	in.read((char*) &initiated, sizeof(bool));
	if(!initiated)
		return in;

	string::size_type nAlphabet, nName;
	DigitalSeq::size_type len;
	string alphabet;

	/* load sizes */
	in.read((char*) &nAlphabet, sizeof(string::size_type));
	in.read((char*) &nName, sizeof(string::size_type));
	in.read((char*) &len, sizeof(DigitalSeq::size_type));

	/* load basic info */
	StringUtils::loadString(alphabet, in, nAlphabet);
	abc = AlphabetFactory::getAlphabetByName(alphabet); /* set alphabet by name */
	StringUtils::loadString(name, in, nName);

	/* load seq */
	StringUtils::loadString(*this, in, len);

	return in;
}

bool DigitalSeq::seqEquals(const string& seq, bool allowDegen) const {
	if(seq.length() != length())
		return false;
	return allowDegen ? seq == toString() : DigitalSeq(abc, name, seq) == *this;
}

} /* namespace EGriceLab */

