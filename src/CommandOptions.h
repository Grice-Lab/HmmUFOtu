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
 * CommandOptions.h
 *
 *  Created on: Jul 15, 2016
 *      Author: zhengqi
 */

#ifndef SRC_COMMANDOPTIONS_H_
#define SRC_COMMANDOPTIONS_H_

#include <string>
#include <set>
#include <map>
#include <vector>

namespace EGriceLab {

using std::string;
using std::map;
using std::vector;

class CommandOptions {
public:
	/** constructors */
	/* construct a CommandOptions from a C main args */
	CommandOptions(int argc, char** argv);

	/** member methods */
	bool hasOpt(const string& name) const {
		return opts.find(name) != opts.end();
	}

	string getOpt(const string& name) const {
		return hasOpt(name) ? opts.find(name)->second : "";
	}

	const char* getOptStr(const string& name) const {
		return hasOpt(name) ? opts.find(name)->second.c_str() : "";
	}

	int numMainOpts() const {
		return mainOpts.size();
	}

	int numOpts() const {
		return opts.size();
	}

	bool empty() const {
		return numMainOpts() == 0 && numOpts() == 0;
	}

	string getMainOpt(int i) const {
		return mainOpts.at(i);
	}

	const string& getProg() const {
		return prog;
	}

	const string& getOptStr() const {
		return optStr;
	}

	string getCmdStr() const {
		return prog + " " + optStr;
	}

private:
	string prog; /* program called */
	string optStr; /* all options as a string */
	vector<string> mainOpts; /* mandatory options not following any -tag-name */
	map<string, string> opts; /* optional named options in -tag-name [value] pairs */
};

} /* namespace EGriceLab */

#endif /* SRC_COMMANDOPTIONS_H_ */
