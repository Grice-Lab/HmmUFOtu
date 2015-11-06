/*
 * DegenAlphabet.cpp
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#include <cassert>
#include <cstdlib>
#include "DegenAlphabet.h"
#include "StringUtils.h"

namespace EGriceLab {
using namespace std;

DegenAlphabet::DegenAlphabet(const string& name, const string& sym_str, const string& synon_str,
			const map<char, string>& my_map, const string& gap) :
				EGriceLab::Alphabet(name, sym_str),/* invoke base class constructor */
				synon(EGriceLab::remove_dup_chars(synon_str)), degen_map(my_map), gap(gap), gap_map() /* zero-initiation */ {
	assert(synon.length() == degen_map.size());
	// init the synon_map
	std::fill_n(synon_map, INT8_MAX + 1, invalid_synon); // fill synon_map with -2
	// set the synon_map
	for(map<char, string>::const_iterator it = degen_map.begin(); it != degen_map.end(); ++it)
		synon_map[it->first] = Alphabet::encode(it->second[::rand() % it->second.length()]); /* set synom map to a random symbol coding */

	// set the gap_map
	for(string::const_iterator it = gap.begin(); it != gap.end(); ++it)
		gap_map[*it] = 1;
}

bool operator==(const DegenAlphabet& lhs, const DegenAlphabet& rhs) {
	if(lhs == rhs) /* same object */
		return true;
	return static_cast<Alphabet> (lhs) == static_cast<Alphabet> (rhs) && lhs.synon == rhs.synon &&
			lhs.degen_map == rhs.degen_map;
}

} /* namespace EGriceLab */

