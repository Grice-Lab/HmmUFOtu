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
 * DegenAlphabet.cpp
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>
#include "DegenAlphabet.h"
#include "StringUtils.h"

namespace EGriceLab {
namespace HmmUFOtu {

using namespace std;

const int8_t DegenAlphabet::INVALID_BASE = -1;
const int8_t DegenAlphabet::GAP_BASE = -2; /* encoded gap symbol */

DegenAlphabet::DegenAlphabet(const string& name, const string& sym_str, const string& synon_str,
			const map<char, string>& my_map, const string& gap) :
				name(name), symbol(StringUtils::remove_dup_chars(sym_str)),
				synon(StringUtils::remove_dup_chars(synon_str)), degen_map(my_map), gap(gap) { /* gapCh default initiated */
	assert(symbol.length() <= INT8_MAX + 1);
	assert(synon.length() == degen_map.size());
	if(!gap.empty())
		gapCh = gap.front();

	// init the sym_map
	std::fill_n(sym_map, INT8_MAX + 1, INVALID_BASE);
	// set the symbol map for both upper and lower cases
	for(int8_t i = 0; i != symbol.length(); ++i) {
		char c = symbol[i];
		assert(std::isupper(c));
		sym_map[c] = i;
		sym_map[std::tolower(c)] = i;
	}

	/* process and update degen_map */
	for(const map<char, string>::value_type& pair : degen_map) { /* set synom map for both upper and lower case symbols */
		char s = pair.first;
		const string& synon = pair.second;
		char c = synon.front(); // use the first synon char
		assert(std::isupper(s) && std::isupper(c));
		/* update degen_map to include lower case */
		degen_map[::tolower(s)] = synon; // lower-case synon still map to upper case symbols
		/* add synon to sym_map */
		sym_map[s] = encode(c);
		sym_map[std::tolower(s)] = encode(c);
	}

	// set the gap_sym
	for(char c : gap)
		sym_map[c] = GAP_BASE;
}

bool DegenAlphabet::isMatch(char c1, char c2) const {
	bool isSynon1 = isSynonymous(c1);
	bool isSynon2 = isSynonymous(c2);
	if(! isSynon1 && ! isSynon2)
		return encode(c1) == encode(c2);
	else if(! isSynon1 && isSynon2)
		return getSynonymous(c2).find(c1) != std::basic_string<int8_t>::npos;
	else if(isSynon1 && ! isSynon2)
		return getSynonymous(c1).find(c2) != std::basic_string<int8_t>::npos;
	else
		return ! StringUtils::common(getSynonymous(c1), getSynonymous(c2)).empty();
}

bool DegenAlphabet::isMatch(char c, int8_t b) const {
	return !isSynonymous(c) && encode(c) == b || /* is not a snynom */
			getSynonymous(c).find(decode(b)) != string::npos; /* search synom */
}

bool operator==(const DegenAlphabet& lhs, const DegenAlphabet& rhs) {
	return lhs.symbol == rhs.symbol && lhs.synon == rhs.synon &&
			lhs.degen_map == rhs.degen_map && lhs.gap == rhs.gap;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

