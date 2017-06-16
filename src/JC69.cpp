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
 * JC69.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: zhengqi
 */

#include "JC69.h"
#include "ProgLog.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const string JC69::name = "JC69";
const Vector4d JC69::pi = Vector4d::Constant(1.0 / 4);

istream& JC69::read(istream& in) {
	string line, tag, value;
	while(in >> tag) {
		if(tag[0] == '#') { /* comment or header */
			std::getline(in, line); /* ignore the entire line */
			continue;
		}
		if(tag == "Type:") {
			in >> value; // read in model type
			if(value != modelType()) {
				errorLog << "Unmatched Model Type!" << endl;
				errorLog << "Trying to read in a " << value << " model into a " << modelType() << " object" << endl;
				in.setstate(ios_base::badbit);
				return in;
			}
			std::getline(in, line); /* ignore the reset of line */
			break;
		}
		else {
			errorLog << "Un-recognized line found in JC69 Model input: tag: " << tag << endl << line << endl;
			in.setstate(ios_base::badbit);
			return in;
		}
	}

	return in;
}

ostream& JC69::write(ostream& out) const {
	out << "# DNA Substitution Model" << endl;
	out << "Type: " << modelType() << endl;

	return out;
}

} /* namespace EGriceLab */
