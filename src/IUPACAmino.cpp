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

#include <string>
#include <iostream>
#include "IUPACAmino.h"

namespace EGriceLab {
namespace HmmUFOtu {

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

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

