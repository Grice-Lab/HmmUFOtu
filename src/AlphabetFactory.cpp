/*
 * AlphabetFactory.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: zhengqi
 */

#include "AlphabetFactory.h"
#include "IUPACNucl.h"
#include "IUPACAmino.h"

#include "StringUtils.h"

namespace EGriceLab {
using namespace std;

const DegenAlphabet* AlphabetFactory::nuclAbc = new IUPACNucl();
const DegenAlphabet* AlphabetFactory::aminoAbc = new IUPACAmino();

const DegenAlphabet* AlphabetFactory::getAlphabetByName(const string& alphabet) {
	string name = StringUtils::toLower(alphabet);
	if(name == "dna" || name == "rna" || alphabet == "IUPACNucl")
		return nuclAbc;
	else if(name == "protein" || alphabet == "IUPACAmino")
		return aminoAbc;
	else
		throw invalid_argument("Unknown alphabet name found: " + alphabet);
}

} /* namespace EGriceLab */

