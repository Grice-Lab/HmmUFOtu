/*
 * AlphabetFactory.h
 *
 *  Created on: Jul 22, 2015
 *      Author: zhengqi
 */

#ifndef SEQCOMMONS_H_
#define SEQCOMMONS_H_
#include <stdexcept>
#include "DegenAlphabet.h"
#include "IUPACNucl.h"
#include "IUPACAmino.h"

namespace EGriceLab {

/**
 * A class for storing common static objects for BioSeq related classes
 */
class AlphabetFactory {
public:
	/* static methods */
	static const DegenAlphabet* getAlphabetByName(const string& alphabet);

private:
	static const DegenAlphabet* nuclAbc;
	static const DegenAlphabet* aminoAbc;
};

} /* namespace EGriceLab */

#endif /* SEQCOMMONS_H_ */
