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
using namespace std;

namespace EGriceLab {

DigitalSeq::DigitalSeq(const DegenAlphabet* abc, const string& name, const string& str) :
				abc(abc), name(name) {
	for(string::const_iterator it = str.begin(); it != str.end(); ++it)
		if(abc->isValid(::toupper(*it)))
			push_back(abc->encode(::toupper(*it))); // use encoded values
}

string DigitalSeq::toString() const {
	string str;
	for(DigitalSeq::const_iterator it = begin(); it != end(); ++it)
		str.push_back(abc->decode(*it));
	return str;
}

DigitalSeq EGriceLab::DigitalSeq::revcom() const {
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
	for(string::const_iterator it = str.begin(); it != str.end(); ++it)
		if(abc->isValid(::toupper(*it)))
			push_back(abc->encode(::toupper(*it)));
	return *this;
}

/*
 * compare two DigitalSeq
 * return true if and only if all residuals are equal and are the same Alphabet
 */
bool operator==(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	if(*(lhs.abc) != *(rhs.abc))
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

DigitalSeq::DigitalSeq(const PrimarySeq& seq) :
	abc(seq.getDegenAlphabet()), name(seq.getId()) {
	append(seq.getSeq());
}

ostream& operator<<(ostream& os, const DigitalSeq& dSeq) {
	for(DigitalSeq::const_iterator it = dSeq.begin(); it != dSeq.end(); ++it)
		os << dSeq.abc->decode(*it);
	return os;
}

} /* namespace EGriceLab */

