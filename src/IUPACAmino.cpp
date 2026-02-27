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
	/* set upper case synonymous */
	IUPAC_map['B'] = string("DN");
	IUPAC_map['X'] = string("ACDEFGHIKLMNPQRSTVWY");
	IUPAC_map['Z'] = string("EQ");
	/* set lower case synonymous */
	IUPAC_map['b'] = string("dn");
	IUPAC_map['x'] = string("acdefghiklmnpqrstvwy");
	IUPAC_map['z'] = string("eq");
	return IUPAC_map;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

