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
 * CommandOptions.cpp
 *
 *  Created on: Jul 15, 2016
 *      Author: zhengqi
 */

#include "CommandOptions.h"
#include <iostream>

namespace EGriceLab {

CommandOptions::CommandOptions(int argc, char** argv) : prog(argv[0])
{
	/* parse options */
	for(int i = 1; i < argc; ++i) {
		optStr += i < argc - 1 ? argv[i] + string(" ") : argv[i];
		if(*argv[i] == '-') { /* a tag name */
			if(i < argc - 1 && *argv[i+1] != '-') {/* a tag value */
				opts[argv[i]] = argv[i+1];
				i++;
			}
			else /* a flag tag */
				opts[argv[i]].push_back('\0'); /* append null values */
		}
		else /* a main opt */
			mainOpts.push_back(argv[i]);
	}
}

} /* namespace EGriceLab */

