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
 * IUPACNucl.h
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#ifndef IUPACNUCL_H_
#define IUPACNUCL_H_

#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include <cctype>
#include "DegenAlphabet.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::string;
using std::map;

class IUPACNucl: public DegenAlphabet {
public:
	/* Constructors */
	/* default constructor */
	IUPACNucl();

	/* destructor, do nothing */
	virtual ~IUPACNucl() { };

	/* member methods */

	/**
	 * Get alias of this alphabet
	 * @override  base class method
	 */
	virtual string getAlias() const {
		return "DNA";
	}

	/**
	 * tell whether has complement symbols, always true
	 * @override  base class method
	 */
	virtual bool hasComplement() const {
		return true;
	}
	/**
	 * Get the complement char of given symbol
	 * @return the complement symbol of matched case, or leave unchanged if not defined
	 */
	virtual char getComplementSymbol(char c) const {
		return !::islower(c) ? compl_map[c] : ::tolower(compl_map[::toupper(c)]);
	}

private:
	/* static initialization method */
	static map<char, string> init_IUPAC_map();
/*	static const map<char, string> degen_map;*/
	int8_t compl_map[INT8_MAX + 1];
};

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* IUPACNUCL_H_ */
