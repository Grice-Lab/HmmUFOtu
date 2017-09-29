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
 * K80.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: zhengqi
 */

#include <iomanip>
#include <cfloat>
#include "K80.h"
#include "ProgLog.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const string K80::name = "K80";
const Vector4d K80::pi = Vector4d::Constant(1.0 / 4);

istream& K80::read(istream& in) {
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
		}
		else if(tag == "kappa:") {
			in >> kappa;
			std::getline(in, line); /* ignore the entire line */
			break;
		}
		else {
			errorLog << "Un-recognized line found in K80 Model input: tag: " << tag << endl << line << endl;
			in.setstate(ios_base::badbit);
			return in;
		}
	}

	setBeta();
	return in;
}

ostream& K80::write(ostream& out) const {
	out << "# DNA Substitution Model" << endl;
	out << "Type: " << modelType() << endl;
	out << std::setprecision(DBL_DIG) << "kappa: " << kappa << endl;

	return out;
}

void K80::trainParams(const vector<Matrix4d>& Pv, const Vector4d& f) {
	/* estimate kappa */
	double Ti = 0, Tv = 0;
	for(vector<Matrix4d>::const_iterator P = Pv.begin(); P != Pv.end(); ++P) {
		Ti += (*P)(A, G) + (*P)(G, A) + (*P)(C, T) + (*P)(T, C);
		Tv += (*P)(A, C) + (*P)(A, T) + (*P)(C, A) + (*P)(C, G) + (*P)(G, C) + (*P)(G, T) + (*P)(T, A) + (*P)(T, G);
	}
	kappa = Ti / Tv;
}

} /* namespace EGriceLab */
