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
#include <cctype>
#include <cassert>
#include "DigitalSeq.h"
#include "StringUtils.h"
#include "AlphabetFactory.h"

using namespace std;

namespace EGriceLab {
namespace HmmUFOtu {

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
	string str(length(), '\0'); // construct a str with enough size
	std::transform(begin(), end(), str.begin(),
	[&](int8_t s) -> char { return abc->decode(s); }); // apply decode transform
	return str;
}

DigitalSeq DigitalSeq::revcom() const {
	if(!abc->hasComplement())
		throw std::invalid_argument("Sequence alphabet " + abc->getName() + " does not support reverse-complement");
	DigitalSeq rcSeq(*this); // make copy of this seq
	std::reverse(rcSeq.begin(), rcSeq.end()); // reverse
	std::transform(rcSeq.begin(), rcSeq.end(), rcSeq.begin(),
			[&] (DigitalSeq::value_type b) -> DigitalSeq::value_type { return abc->encode(abc->getComplementSymbol(abc->decode(b))); }
	); // apply complement transform
	return rcSeq;
}

string DigitalSeq::join(const string& sep) {
	string str;
	str.reserve(2 * length()); // make enough reserve
	for(DigitalSeq::value_type b : *this) {
		if(!str.empty())
			str += sep;
		str += abc->decode(b);
	}
	return str;
}

DigitalSeq& DigitalSeq::append(const string& str) {
	for(string::value_type c : str) {
		c = std::toupper(c);
		if(abc->isValid(c))
			push_back(abc->encode(c));
	}
	return *this;
}

ostream& DigitalSeq::save(ostream& out, bool withAbc) const {
	/* save flag */
	bool flag = abc != NULL && withAbc;
	out.write((const char*) &flag, sizeof(bool));

	/* save alphabet, if requested */
	if(flag)
		StringUtils::saveString(abc->getName(), out);

	/* save data */
	StringUtils::saveString(name, out);
	StringUtils::saveString(*this, out);
	return out;
}

istream& DigitalSeq::load(istream& in) {
	/* load flag */
	bool flag;
	in.read((char*) &flag, sizeof(bool));

	/* load alphabet, if requested */
	if(flag) {
		string alphabet;
		StringUtils::loadString(alphabet, in);
		abc = AlphabetFactory::getAlphabetByName(alphabet);
	}

	/* load data */
	StringUtils::loadString(name, in);
	StringUtils::loadString(*this, in);
	return in;
}

bool DigitalSeq::seqEquals(const string& seq, bool allowDegen) const {
	if(seq.length() != length())
		return false;
	return allowDegen ? seq == toString() : DigitalSeq(abc, name, seq) == *this;
}

ostream& operator<<(ostream& os, const DigitalSeq& dSeq) {
	for(DigitalSeq::value_type b : dSeq)
		os << dSeq.abc->decode(b);
	return os;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
