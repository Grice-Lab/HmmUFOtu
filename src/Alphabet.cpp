/*
 * Alphabet.cpp
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */


#include "Alphabet.h"
#include "StringUtils.h"
#include <algorithm>
#include <cassert>
using namespace std;

namespace EGriceLab {

Alphabet::Alphabet(const string& name, const string& sym_str) :
		name(name), sym(EGriceLab::remove_dup_chars(sym_str)) {
	assert(sym.length() <= INT8_MAX + 1);
	// init the inmap
	std::fill_n(sym_map, INT8_MAX + 1, invalid_sym); // fill inmap with -1
	// set the map
	for(int8_t i = 0; i != sym.length(); ++i)
		sym_map[sym[i]] = i;
}

} /* namespace EGriceLab */

