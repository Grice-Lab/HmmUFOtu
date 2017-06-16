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
 * IUPACAmino.h
 *
 *  Created on: May 5, 2015
 *      Author: zhengqi
 */

#ifndef IUPACAMINO_H_
#define IUPACAMINO_H_

#include <map>
#include <string>
#include <stdexcept>
#include "DegenAlphabet.h"

namespace EGriceLab {
using std::string;
using std::map;

class IUPACAmino: public DegenAlphabet {
public:
	/* Constructors */
	/* default constructor */
	IUPACAmino() : DegenAlphabet("IUPACAmino",
			"ACDEFGHIKLMNPQRSTVWY", "BXZ", init_IUPAC_map()) {
	}

	/* destructor, do nothing */
	virtual ~IUPACAmino() { };

	/* member methods */
	/**
	 * Get alias of this alphabet
	 * @override  base class method
	 */
	virtual string getAlias() const {
		return "AMINO";
	}

	/**
	 * check whether has complement, always false
	 * @override  base class method
	 */
	virtual bool hasComplement() const {
		return false;
	}

	/**
	 * Get the complement char of given symbol
	 * @return unchanged amino acids don't have complementary symbols
	 */
	virtual char getComplementSymbol(char c) const {
		return c;
	}

private:
	/* static initialization method */
	static map<char, string> init_IUPAC_map();
};

} /* namespace EGriceLab */

#endif /* IUPACAMINO_H_ */
