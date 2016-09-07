/*
 * PhyloTreeNode.cpp
 *
 *  Created on: Mar 25, 2016
 *      Author: zhengqi
 */

#include <set>
#include <stack>
#include <iostream>
#include <fstream>
#include "PhyloTree.h"

namespace EGriceLab {
using namespace std;

int PhyloTree::numNodes() const {
	int n = 0;
	set<const PT*> visited;
	stack<const PT*> S;

	S.push(this);
	while(!S.empty()) {
		const PT* v = S.pop();
		if(visited.find(v) != visited.end()) {
			visited.insert(v);
			for(vector<PT>::const_iterator it = children.begin(); it != children.end(); ++it) {
				const PT* p = &*it;
				S.push(p);
			}
		}
	}
	return visited.size();
}

int PhyloTree::readNiwickTree(const string& treefn) {
	namespace qi = boost::spirit::qi;

	ifstream in(treefn);
	if(!in.is_open()) {
		cerr << "Error: Could not open input file: " << treefn << endl;
		return -1;
	}

	EGriceLab::newick_grammar<istream_iterator<char> > grammar;
	EGriceLab::PhyloTree tree;

	std::istream_iterator<char> iter(in);
	std::istream_iterator<char> end;
	// parse
	bool result = qi::phrase_parse(iter, end, grammar, qi::space, tree);

	if(result && iter == end)
		return numNodes();
	else {
		cout << "Parsing newick tree failed." << endl;
		return -1;
	}
}

} /* namespace EGriceLab */

