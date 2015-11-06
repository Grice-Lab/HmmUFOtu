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
 * compare two DigitalSeq strict weak order, at lexical order
 * return true if and only if lhs is strictly less than rhs
 */
bool operator<(const DigitalSeq& lhs, const DigitalSeq& rhs) {
	// construct two temporary C-strs
	char lstr[lhs.length() + 1];
	char rstr[rhs.length() + 1];
	for(DigitalSeq::size_type i = 0; i != lhs.length(); ++i)
		lstr[i] = static_cast<char>(lhs[i]);
	for(DigitalSeq::size_type i = 0; i != rhs.length(); ++i)
		rstr[i] = static_cast<char>(rhs[i]);
	lstr[lhs.length()] = '\0'; // null terminal
	rstr[rhs.length()] = '\0'; // null terminal
	return string(lstr) < string(rstr);
}

ostream& operator<<(ostream& os, const DigitalSeq& dSeq) {
	for(DigitalSeq::const_iterator it = dSeq.begin(); it != dSeq.end(); ++it)
		os << dSeq.abc->decode(*it);
	return os;
}

} /* namespace EGriceLab */

EGriceLab::DigitalSeq::DigitalSeq(const DegenAlphabet& dgAbc, const string& name, const string& str) :
		abc(&dgAbc), name(name) {
	for(string::const_iterator it = str.begin(); it != str.end(); ++it) {
		int ds = abc->encode(toupper(*it));
		if(ds >= 0)
			push_back(ds);
	}
}

string EGriceLab::DigitalSeq::toString() const {
	string str;
	for(DigitalSeq::const_iterator it = begin(); it != end(); ++it)
		str.push_back(abc->decode(*it));
	return str;
}

EGriceLab::DigitalSeq EGriceLab::DigitalSeq::revcom() const {
	if(!abc->hasComplement())
		throw std::invalid_argument("Sequence alphabet " + abc->getName() + " does not support reverse-complement");
	DigitalSeq revSeq(*abc, name); // make an empty copy with same DegebAlphabet and name
	for(DigitalSeq::const_iterator it = begin(); it != end(); ++it)
		revSeq.push_back(abc->encode(abc->getComplementSymbol(abc->decode(*it))));
	return revSeq;
}

string EGriceLab::DigitalSeq::join(const string& sep) {
	ostringstream ostr;
	for(const_iterator it = begin(); it != end(); ++it) {
		if(it != begin())
			ostr << sep;
		ostr << *it;
	}
	return ostr.str();
}

EGriceLab::DigitalSeq& EGriceLab::DigitalSeq::append(const string& str) {
	for(string::const_iterator it = str.begin(); it != str.end(); ++it) {
		int ds = abc->encode(toupper(*it));
		if(ds != string::npos)
			push_back(ds);
	}
	return *this;
}


istream& EGriceLab::DigitalSeq::readFASTANext(DigitalSeq& ds, istream& is) {
	string line;
	ds.clear(); // clear old seq
	while(is.peek() != '>') // jump to the next record
		getline(is, line);

	while((ds.empty() || is.peek() != '>') && getline(is, line)) {
		if(line.empty())
			continue; // ignore empty lines
		if(line[0] == '>') { // def line
			istringstream iss(line.substr(1)); // ignore the leading '>'
			string name;
			iss >> name;
			ds.setName(name);
		}
		else {
			ds.append(line);
		}
	}
	return is;
}

istream& EGriceLab::DigitalSeq::readFASTQNext(DigitalSeq& ds, istream& is) {
	string line;
	string name;
	ds.clear();
	if(is.peek() != '@') // jump to the next record
		getline(is, line);

	/* read def line */
	getline(is, line);
	istringstream iss(line.substr(1)); // ignore the leading '@'
	iss >> name;
	ds.setName(name);

	/* read seq line */
	getline(is, line);
	ds.append(line);

	/* Ignore sep line and qual line */
	getline(is, line);
	getline(is, line);
	return is;
}
