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
#include "DegenAlphabet.h"
#include "StringUtils.h"
#include <iostream>

namespace EGriceLab {
namespace HmmUFOtu {

using namespace std;

const int8_t DegenAlphabet::INVALID_BASE = -1;
const int8_t DegenAlphabet::GAP_BASE = -2; /* encoded gap symbol */

DegenAlphabet::DegenAlphabet(const string& name, const string& sym_str, const string& synon_str,
			const map<char, string>& my_map, const string& gap) :
				name(name), symbol(StringUtils::remove_dup_chars(sym_str)),
				synon(StringUtils::remove_dup_chars(synon_str)), degen_map(my_map), gap(gap) {
	assert(symbol.length() <= INT8_MAX + 1);
	assert(synon.length() == degen_map.size());
	!gap.empty() ? gapCh = gap[0] : DEFAULT_GAP_CHAR;

	// init the sym_map
	std::fill_n(sym_map, INT8_MAX + 1, INVALID_BASE);
	// set the symbol map
	for(int8_t i = 0; i != symbol.length(); ++i)
		sym_map[symbol[i]] = i;

	// set the synon_map
	for(map<char, string>::const_iterator it = degen_map.begin(); it != degen_map.end(); ++it)
		sym_map[it->first] = encode(it->second[0]); /* set synom map to the first symbol */

	// set the gap_sym
	for(string::const_iterator it = gap.begin(); it != gap.end(); ++it)
		sym_map[*it] = GAP_BASE;
}

bool DegenAlphabet::isMatch(char c1, char c2) const {
	return StringUtils::common(c1 + getSynonymous(c1), c2 + getSynonymous(c2)) > 0;
}

bool DegenAlphabet::isMatch(char c, int8_t b) const {
	string synon = c + getSynonymous(c);
	for(string::const_iterator ch = synon.begin(); ch != synon.end(); ++ch)
		if(encode(*ch) == b)
			return true;
	return false;
}

bool operator==(const DegenAlphabet& lhs, const DegenAlphabet& rhs) {
	return lhs.symbol == rhs.symbol && lhs.synon == rhs.synon &&
			lhs.degen_map == rhs.degen_map && lhs.gap == rhs.gap;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

