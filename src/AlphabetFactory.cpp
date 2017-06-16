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

