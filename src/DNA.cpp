/*
 * DNA.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: zhengqi
 */

#include "DNA.h"

namespace EGriceLab {

map<char, string> DNA::init_DNA_map() {
	map<char, string> dna_map;
	dna_map['U'] = string("T");
	dna_map['N'] = string("ACGT");
	return dna_map;
}

DNA::DNA() : EGriceLab::DegenAlphabet("DNA", "ACGT", "UN", init_DNA_map()), compl_map() /* zero-initiazation */ {
	// std::cerr << "Constructing IUPACNucl" << std::endl;
	compl_map['A'] = 'T';
	compl_map['T'] = 'A';
	compl_map['C'] = 'G';
	compl_map['G'] = 'C';
	compl_map['U'] = 'A';
	compl_map['N'] = 'N';
}

} /* namespace EGriceLab */
