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
 * DNASubModelFactory.h
 *  A factory class with static methods creating new empty DNASubstitution models base on type names
 *  Created on: Dec 16, 2016
 *      Author: zhengqi
 */

#ifndef SRC_DNASUBMODELFACTORY_H_
#define SRC_DNASUBMODELFACTORY_H_

#include <stdexcept>
#include <string>
#include "DNASubModel.h"


namespace EGriceLab {

using std::string;

class DNASubModelFactory {
public:
	static DNASubModel* createModel(const string& type);
};

} /* namespace EGriceLab */

#endif /* SRC_DNASUBMODELFACTORY_H_ */
