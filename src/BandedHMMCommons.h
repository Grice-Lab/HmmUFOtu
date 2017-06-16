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
 * BandedHMMCommons.h
 *	Contains basic public class for BandedHMM algorithm
 *  Created on: Aug 3, 2015
 *      Author: zhengqi
 */

#ifndef BANDEDHMMCOMMONS_H_
#define BANDEDHMMCOMMONS_H_
#include <string>

namespace EGriceLab {
using std::string;
/**
 * A public class for describing a region on the concensus seq (CS)
 */
struct CSLoc {
	/* constructors */
	/**
	 * Default constructor
	 */
	CSLoc() : start(0), end(0), prob(0) { }

	/**
	 * Construct a CSLoc at given loc
	 */
	CSLoc(int start, int end, const string& CS, int prob = 0)
		: start(start), end(end), CS(CS), prob(prob) { }

	int start;
	int end;
	string CS;
	int prob;
};

}

#endif /* BANDEDHMMCOMMONS_H_ */
