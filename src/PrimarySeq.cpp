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
		char rCh = abc->getComplementSymbol(::toupper(*it));
		revcomSeq.push_back(::isupper(*it) ? rCh : ::tolower(rCh));
	}
	return PrimarySeq(abc, id, revcomSeq, desc);
}

string::size_type PrimarySeq::numGap() const {
	string::size_type n = 0;
	for(string::const_iterator it = seq.begin(); it != seq.end(); ++it)
		if(abc->isGap(*it))
			n++;
	return n;
}

} /* namespace EGriceLab */
