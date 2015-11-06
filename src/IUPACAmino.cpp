/*
 * DegenAlphabet.cpp
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include "IUPACAmino.h"

namespace EGriceLab {
map<char, string> IUPACAmino::init_IUPAC_map() {
	map<char, string> IUPAC_map;
	//IUPAC_map['A'] = string("A");
	IUPAC_map['B'] = string("DN");
	//IUPAC_map['C'] = string("C");
	//IUPAC_map['D'] = string("D");
	//IUPAC_map['E'] = string("E");
	//IUPAC_map['F'] = string("F");
	//IUPAC_map['G'] = string("G");
	//IUPAC_map['H'] = string("H");
	//IUPAC_map['I'] = string("I");
	//IUPAC_map['K'] = string("K");
	//IUPAC_map['L'] = string("L");
	//IUPAC_map['M'] = string("M");
	//IUPAC_map['N'] = string("N");
	//IUPAC_map['P'] = string("P");
	//IUPAC_map['Q'] = string("Q");
	//IUPAC_map['R'] = string("R");
	//IUPAC_map['S'] = string("S");
	//IUPAC_map['T'] = string("T");
	//IUPAC_map['V'] = string("V");
	//IUPAC_map['W'] = string("W");
	IUPAC_map['X'] = string("ACDEFGHIKLMNPQRSTVWY");
	//IUPAC_map['Y'] = string("Y");
	IUPAC_map['Z'] = string("EQ");
	return IUPAC_map;
}

} /* namespace EGriceLab */

