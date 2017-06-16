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
 * DNASubModelFactory.cpp
 *
 *  Created on: Dec 16, 2016
 *      Author: zhengqi
 */

#include "DNASubModelFactory.h"
#include "GTR.h"
#include "TN93.h"
#include "HKY85.h"
#include "F81.h"
#include "K80.h"
#include "JC69.h"

namespace EGriceLab {

DNASubModel* DNASubModelFactory::createModel(const string& type) {
	if(type == "GTR")
		return new GTR();
	else if(type == "TN93")
		return new TN93();
	else if(type == "HKY85")
		return new HKY85();
	else if(type == "F81")
		return new F81();
	else if(type == "K80")
		return new K80();
	else if(type == "JC69")
		return new JC69();
	else
		throw std::invalid_argument("Unknown DNA Substitution Model type: " + type);
}

} /* namespace EGriceLab */


