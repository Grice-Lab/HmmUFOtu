/*
 * SeqCommons.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: zhengqi
 */

#include "SeqCommons.h"
#include "StringUtils.h"

namespace EGriceLab {
using namespace std;

const IUPACNucl* const SeqCommons::nuclAbc = new IUPACNucl();
const IUPACAmino* const SeqCommons::aminoAbc = new IUPACAmino();

const DegenAlphabet* SeqCommons::getAlphabetByName(const string& alphabet) {
	const string& name = toLower(alphabet);
	if(name == "dna" || name == "rna" || alphabet == "IUPACNucl")
		return nuclAbc;
	else if(name == "protein" || alphabet == "IUPACAmino")
		return aminoAbc;
	else
		throw invalid_argument("Unknown alphabet name found: " + alphabet);
}

} /* namespace EGriceLab */

