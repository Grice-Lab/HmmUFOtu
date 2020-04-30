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
 * NewickTree.cpp
 *
 *  Created on: Dec 2, 2016
 *      Author: zhengqi
 */

#include "StringUtils.h"
#include "NewickTree.h"

namespace EGriceLab {
namespace HmmUFOtu {

using namespace std;

const string& NewickTree::INVALID_CHARS = "()[]':;,";

istream& NewickTree::read(istream& in) {
	namespace qi = boost::spirit::qi;

	string content; /* store entire newick file in a string */

	newick_grammar<std::string::const_iterator> grammar;

	/* copy the whole content in input */
	in.unsetf(std::ios::skipws); /* do not skip whitespaces */
	std::copy(std::istream_iterator<char>(in), std::istream_iterator<char>(), std::back_inserter(content));

	string::const_iterator iter = content.begin();
	string::const_iterator end = content.end();
	// clear old data
	clear();
	// parse
	bool result = qi::phrase_parse(iter, end, grammar, qi::space, *this);

	if(!(result && iter == end))
		in.setstate(ios_base::badbit);

	return in;
}

ostream& NewickTree::write(ostream& out) const {
	bool first = true;
	if(!children.empty()) {
		out << '(';
		for(std::vector<NT>::const_iterator it = children.begin(); it != children.end(); ++it) {
			out << (first ? "" : ",");
			it->write(out);
			first = false;
		}
		out << ')';
	}
	out << quoteName(name);
	if(length >= 0)
		out << ':' << length;

	return out;
}

bool NewickTree::isNewickFileExt(const string& fn) {
	return StringUtils::endsWith(fn, ".tree") || StringUtils::endsWith(fn, ".tre");
}


} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
