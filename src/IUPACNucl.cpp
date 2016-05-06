/*
 * DegenAlphabet.cpp
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include "IUPACNucl.h"

namespace EGriceLab {

map<char, string> IUPACNucl::init_IUPAC_map() {
	map<char, string> IUPAC_map;
	IUPAC_map['U'] = string("T");
	IUPAC_map['M'] = string("AC");
	IUPAC_map['R'] = string("AG");
	IUPAC_map['W'] = string("AT");
	IUPAC_map['S'] = string("CG");
	IUPAC_map['Y'] = string("CT");
	IUPAC_map['K'] = string("GT");
	IUPAC_map['V'] = string("ACG");
	IUPAC_map['H'] = string("ACT");
	IUPAC_map['D'] = string("AGT");
	IUPAC_map['B'] = string("CGT");
	IUPAC_map['N'] = string("ACGT");
	return IUPAC_map;
}

IUPACNucl::IUPACNucl() : EGriceLab::DegenAlphabet("IUPACNucl", "ACGT", "UMRWSYKVHDBN", init_IUPAC_map()) {
	/* init compl_map with self complementary */
	for(int8_t i = 0; i <= INT8_MAX; ++i)
		compl_map[i] = i;
	/* override DNA complements */
	compl_map['A'] = 'T';
	compl_map['T'] = 'A';
	compl_map['C'] = 'G';
	compl_map['G'] = 'C';
	compl_map['U'] = 'A';
	compl_map['Y'] = 'R';
	compl_map['R'] = 'Y';
	compl_map['S'] = 'S';
	compl_map['W'] = 'W';
	compl_map['K'] = 'M';
	compl_map['M'] = 'K';
	compl_map['B'] = 'V';
	compl_map['V'] = 'B';
	compl_map['D'] = 'H';
	compl_map['H'] = 'D';
	compl_map['N'] = 'N';
}

} /* namespace EGriceLab */

