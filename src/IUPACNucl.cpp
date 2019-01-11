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
#include "IUPACNucl.h"

namespace EGriceLab {
namespace HmmUFOtu {

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

IUPACNucl::IUPACNucl() : DegenAlphabet("IUPACNucl", "ACGT", "UMRWSYKVHDBN", init_IUPAC_map()) {
	/* init compl_map with self complementary */
	for(int8_t i = 0; i != INT8_MAX; ++i)
		compl_map[i] = i;
	/* upper case complements */
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

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

