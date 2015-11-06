/*
 * SeqCommons.h
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
class SeqCommons {
public:
	/* shared static objects intented to be used by other class */
	//static const DNA* const dnaAbc;
	static const IUPACNucl* const nuclAbc;
	static const IUPACAmino* const aminoAbc;

	/* static methods */
	static const DegenAlphabet* getAlphabetByName(const string& alphabet);
};

} /* namespace EGriceLab */

#endif /* SEQCOMMONS_H_ */
