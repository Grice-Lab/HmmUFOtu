/*
 * PT_IO_test1.cpp
 *
 *  Created on: Aug 25, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include "PhyloTree.h"
using namespace std;

int main(int argc, const char* argv[]) {

	cerr << "Main entered" << endl;

	namespace qi = boost::spirit::qi;
	std::string str;

	EGriceLab::newick_grammar<std::string::const_iterator> grammar;

	cerr << "grammar initiated" << endl;

	EGriceLab::PT tree;

	cerr << "tree initiated" << endl;

	while(std::getline(std::cin, str)) {
		if(str.empty() || str[0] == 'q' || str[0] == 'Q')
			break;
		// parse
		std::string::const_iterator iter = str.begin();
		std::string::const_iterator end = str.end();
		bool result = qi::phrase_parse(iter, end, grammar, qi::space, tree);

		if(result && iter == str.end()) {
			std::cout << "Parsing succeeded. Tree:" << endl <<
					tree << std::endl;
		}
		else {
			std::cout << "Parsing failed. Stopped at: " << std::string(iter, end) << std::endl;
			return 1;
		}
	}

	std::cout << "Bye... :-)\n\n";
	return 0;
}
