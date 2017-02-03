/*
 * NewickTree.cpp
 *
 *  Created on: Dec 2, 2016
 *      Author: zhengqi
 */

#include "StringUtils.h"
#include "NewickTree.h"

namespace EGriceLab {
using namespace std;

const string& NewickTree::INVALID_CHARS = "()[]':;,";

istream& NewickTree::read(istream& in) {
	namespace qi = boost::spirit::qi;

	string content; /* store entire newick file in a string */

	EGriceLab::newick_grammar<std::string::const_iterator> grammar;

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
	if(StringUtils::containsWhiteSpace(name) || StringUtils::containsAny(name, INVALID_CHARS)) // name contains INVALID CHARS
		out << "'" << name << "'";
	else
		out << name;
	if(length > 0)
		out << ':' << length;

	return out;
}


} /* namespace EGriceLab */
