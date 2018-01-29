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
 * BandedHMMP7Prior.cpp
 *
 *  Created on: Jun 13, 2016
 *      Author: zhengqi
 */

#include "BandedHMMP7Prior.h"

namespace EGriceLab {
namespace HmmUFOtu {

const double BandedHMMP7Prior::DEFAULT_ABS_EPS_COST = 0;
const double BandedHMMP7Prior::DEFAULT_REL_EPS_COST = 1e-6;
const double BandedHMMP7Prior::DEFAULT_ABS_EPS_PARAMS = 0;
const double BandedHMMP7Prior::DEFAULT_REL_EPS_PARAMS = 1e-4;


istream& operator>>(istream& in, BandedHMMP7Prior& pri) {
	string head;
	while(std::getline(in, head)) {
		if(StringUtils::startsWith(head, "Match emission:"))
			in >> pri.dmME;
		else if(StringUtils::startsWith(head, "Insert emission:"))
			in >> pri.dmIE;
		else if(StringUtils::startsWith(head, "Match transition:"))
			in >> pri.dmMT;
		else if(StringUtils::startsWith(head, "Insert transition:"))
			in >> pri.dmIT;
		else if(StringUtils::startsWith(head, "Delete transition:"))
			in >> pri.dmDT;
		else
			continue;
	}
	if(!(pri.dmME.getK() > 0 && pri.dmIE.getK() > 0
			&& pri.dmMT.getK() > 0 && pri.dmIT.getK() > 0 && pri.dmDT.getK() > 0)) {
		std::cerr << "Empty or partial BandedHMMP7Prior input" << endl;
		in.setstate(std::ios_base::failbit);
	}
	return in;
}

ostream& operator<<(ostream& out, const BandedHMMP7Prior& pri) {
	out << "Match emission:" << endl << pri.dmME;
	out << "Insert emission:" << endl << pri.dmIE;
	out << "Match transition:" << endl << pri.dmMT;
	out << "Insert transition:" << endl << pri.dmIT;
	out << "Delete transition:" << endl << pri.dmDT;
	return out;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
