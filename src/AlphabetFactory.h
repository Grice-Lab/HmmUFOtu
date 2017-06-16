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

	static const DegenAlphabet* nuclAbc;
	static const DegenAlphabet* aminoAbc;
};

} /* namespace EGriceLab */

#endif /* SEQCOMMONS_H_ */
