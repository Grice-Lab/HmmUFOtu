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
using namespace std;

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
		sym_map[it->first] = encode(it->second.front()); /* set synom map to the first symbol */

	// set the gap_sym
	for(string::const_iterator it = gap.begin(); it != gap.end(); ++it)
		sym_map[*it] = GAP_BASE;
}

bool operator==(const DegenAlphabet& lhs, const DegenAlphabet& rhs) {
	return lhs.symbol == rhs.symbol && lhs.synon == rhs.synon &&
			lhs.degen_map == rhs.degen_map && lhs.gap == rhs.gap;
}

} /* namespace EGriceLab */

