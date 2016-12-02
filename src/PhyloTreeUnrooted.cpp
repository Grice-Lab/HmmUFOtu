/*
 * PhyloTreeUnrooted.cpp
 *
 *  Created on: Dec 1, 2016
 *      Author: zhengqi
 */


#include "PhyloTreeUnrooted.h"

namespace EGriceLab {
using namespace std;
using namespace EGriceLab;

const double PhyloTreeUnrooted::MIN_EXPONENT = std::numeric_limits<float>::min_exponent;

int PhyloTreeUnrooted::numNodes() const {
	return 0;
}

int PhyloTreeUnrooted::readNiwickTree(const string& treefn) {
	namespace qi = boost::spirit::qi;

	ifstream in(treefn.c_str());
	if(!in.is_open()) {
		cerr << "Error: Could not open input file: " << treefn << endl;
		return -1;
	}

	EGriceLab::newick_grammar<std::string::const_iterator> grammar;

	string content;

	/* copy the whole content in input */
	in.unsetf(std::ios::skipws);
	std::copy(std::istream_iterator<char>(in), std::istream_iterator<char>(), std::back_inserter(content));

	string::const_iterator iter = content.begin();
	string::const_iterator end = content.end();
	// parse
	bool result = qi::phrase_parse(iter, end, grammar, qi::space, *this);

	if(result && iter == end)
		return numNodes();
	else {
		cerr << "Parsing Newick tree failed" << endl;
		return -1;
	}
}

} /* namespace EGriceLab */
