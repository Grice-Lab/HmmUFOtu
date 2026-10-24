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
 * DNA.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: zhengqi
 */

#include "DNA.h"

namespace EGriceLab {
namespace HmmUFOtu {

map<char, string> DNA::init_DNA_map() {
	map<char, string> dna_map;
	// upper case syn_map
	dna_map['U'] = string("T");
	dna_map['N'] = string("ACGT");
	// lower case syn_map
	dna_map['u'] = string("t");
	dna_map['n'] = string("acgt");
	return dna_map;
}

DNA::DNA() : DegenAlphabet("DNA", "ACGT", "UN", init_DNA_map()), compl_map() /* zero-initiazation */ {
	// std::cerr << "Constructing IUPACNucl" << std::endl;
	compl_map['A'] = 'T';
	compl_map['T'] = 'A';
	compl_map['C'] = 'G';
	compl_map['G'] = 'C';
	compl_map['U'] = 'A';
	compl_map['N'] = 'N';
    /* lower case complements */
    compl_map['a'] = 't';
    compl_map['t'] = 'a';
    compl_map['c'] = 'g';
    compl_map['g'] = 'c';
    compl_map['u'] = 'a';
    compl_map['n'] = 'n';
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
