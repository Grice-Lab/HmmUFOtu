/*
 * DigitalSeq.cpp
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include "DigitalSeq.h"

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
	abc(seq.getDegenAlphabet()), name(seq.getId()) {
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

	string alphabet = abc->getName();
	string::size_type nAlphabet = alphabet.length();
	string::size_type nName = name.length();
	DigitalSeq::size_type len = length();

	/* save basic info */
	out.write((const char*) &nAlphabet, sizeof(string::size_type));
	out.write(alphabet.c_str(), nAlphabet + 1); /* write null terminal */
	out.write((const char*) &nName, sizeof(string::size_type));
	out.write(name.c_str(), nName + 1); /* write null terminal */
	out.write((const char*) &len, sizeof(DigitalSeq::size_type));

	/* save encoded seq */
	out.write((const char*) c_str(), (len + 1) * sizeof(int8_t)); /* write terminal null */
	return out;
}

istream& DigitalSeq::load(istream& in) {
	/* read flag */
	bool initiated;
	in.read((char*) &initiated, sizeof(bool));
	if(!initiated)
		return in;

	string::size_type nAlphabet, nName;
	string alphabet;
	DigitalSeq::size_type len;
	char* buf = NULL;

	/* load basic info */
	in.read((char*) &nAlphabet, sizeof(string::size_type));
	buf = new char[nAlphabet + 1];
	in.read(buf, nAlphabet + 1);
	alphabet.assign(buf, nAlphabet);
	delete[] buf;
	abc = AlphabetFactory::getAlphabetByName(alphabet); // set alphabet by name

	in.read((char*) &nName, sizeof(string::size_type));
	buf = new char[nName + 1];
	in.read(buf, nName + 1);
	name.assign(buf, nName);
	delete[] buf;

	in.read((char*) &len, sizeof(DigitalSeq::size_type));

	/* load encoded seq */
	buf = new char[(len + 1) * sizeof(int8_t)];
	in.read(buf, len + 1);
	assign((int8_t*) buf, len); /* override old data */
	delete[] buf;

	return in;
}

} /* namespace EGriceLab */

