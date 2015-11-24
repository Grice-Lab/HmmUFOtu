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

namespace EGriceLab {
using namespace std;

DegenAlphabet::DegenAlphabet(const string& name, const string& sym_str, const string& synon_str,
			const map<char, string>& my_map, const string& gap) :
				name(name), symbol(EGriceLab::remove_dup_chars(sym_str)),
				synon(EGriceLab::remove_dup_chars(synon_str)), degen_map(my_map), gap(gap), gap_map() /* zero-initiation */ {
	assert(symbol.length() <= INT8_MAX + 1);
	assert(synon.length() == degen_map.size());

	// init the sym_map
	std::fill_n(sym_map, INT8_MAX + 1, invalid_sym);
	// set the symbol map
	for(int8_t i = 0; i != symbol.length(); ++i)
		sym_map[symbol[i]] = i;

	// set the synon_map
	for(map<char, string>::const_iterator it = degen_map.begin(); it != degen_map.end(); ++it)
		sym_map[it->first] = encode(it->second[::rand() % it->second.length()]); /* set synom map to a random symbol coding */

	// set the gap_map
	for(string::const_iterator it = gap.begin(); it != gap.end(); ++it)
		gap_map[*it] = 1;
}

bool operator==(const DegenAlphabet& lhs, const DegenAlphabet& rhs) {
	return lhs.symbol == rhs.symbol && lhs.synon == rhs.synon &&
			lhs.degen_map == rhs.degen_map && lhs.gap == rhs.gap;
}

} /* namespace EGriceLab */

