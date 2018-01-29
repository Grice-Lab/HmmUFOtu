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
 * DNA.h
 *
 *  Created on: Oct 27, 2015
 *      Author: zhengqi
 */

#ifndef DNA_H_
#define DNA_H_

#include "DegenAlphabet.h"

namespace EGriceLab {
namespace HmmUFOtu {

class DNA: public DegenAlphabet {
public:
	/* Constructors */
	/* default constructor */
	DNA();

	/* destructor, do nothing */
	virtual ~DNA() { };

	/* member methods */
	/* implementation of abstract superclass methods */
	/**
	 * always return true
	 */
	bool hasComplement() const {
		return true;
	}
	/**
	 * Get the complement char of given symbol
	 * @return the complement symbol, or '\0' if not a valid symbol
	 */
	char getComplementSymbol(char c) const {
		return compl_map[c];
	}

private:
	/* static initialization method */
	static map<char, string> init_DNA_map();
/*	static const map<char, string> degen_map;*/
	//map<char, char> compl_map;
	char compl_map[INT8_MAX + 1];
};


} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* DNA_H_ */
