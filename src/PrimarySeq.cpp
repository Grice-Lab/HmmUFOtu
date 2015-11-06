/*
 * PrimarySeq.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: zhengqi
 */

#include <algorithm>
#include <cctype>
#include "PrimarySeq.h"

namespace EGriceLab {
using namespace std;

bool PrimarySeq::isValidate() const {
	for(string::const_iterator it = seq.begin(); it != seq.end(); ++it)
		if(!abc->isValid(::toupper(*it))) // test synonymous case insensitive
			return false;
	return true;
}

PrimarySeq& PrimarySeq::removeGaps() {
	if(seq.length() == 0)
		return *this;
	// remove gaps backwards
	for(string::size_type i = seq.length(); i != 0; --i) {
		if(abc->isGap(seq[i-1])) {
			seq.erase(i-1, 1);
			if(!qual.empty())
				qual.erase(i-1, 1);
		}
	}
	return *this;
}

PrimarySeq PrimarySeq::revcom() const {
	if(!abc->hasComplement())
		throw logic_error("This seq's alphabet " + abc->getName() + " doesn't support reverse-complement action");
	string revcomSeq;
	for(string::const_reverse_iterator it = seq.rbegin(); it != seq.rend(); ++it) {
		char rCh;
		if(::isupper(*it))
			rCh = abc->getComplementSymbol(*it); // keep in upper-case
		else
			rCh = ::tolower(abc->getComplementSymbol(::toupper(*it))); // keep in lower-case
		revcomSeq.push_back(rCh);
	}
	return PrimarySeq(abc->getName(), id, seq, desc);
}

PrimarySeq PrimarySeq::trunc(string::size_type pos,
		string::size_type len) const {
	PrimarySeq newSeq(*this); // make a copy
	newSeq.setSeq(seq.substr(pos, len));
	if(!qual.empty())
		newSeq.setQual(qual.substr(pos, len));
	return newSeq;
}

string::size_type PrimarySeq::numGap() const {
	string::size_type n = 0;
	for(string::const_iterator it = seq.begin(); it != seq.end(); ++it)
		if(abc->isGap(*it))
			n++;
	return n;
}

} /* namespace EGriceLab */
