/*
 * PhyloTreeNode.cpp
 *
 *  Created on: Mar 25, 2016
 *      Author: zhengqi
 */

#include <set>
#include <stack>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include "PhyloTree.h"

namespace EGriceLab {
using namespace std;

int PhyloTree::numNodes() const {
	int n = 0;
	set<const PT*> visited;
	stack<const PT*> S;

	S.push(this);
	while(!S.empty()) {
		const PT* v = S.top();
		S.pop();
		if(visited.find(v) != visited.end()) {
			visited.insert(v);
			for(vector<PT>::const_iterator it = v->children.begin(); it != v->children.end(); ++it) {
				const PT* p = &*it;
				S.push(p);
			}
		}
	}
	return visited.size();
}

int PhyloTree::numLeaves() const {
	int n = 0;
	set<const PT*> visited;
	stack<const PT*> S;

	S.push(this);
	while(!S.empty()) {
		const PT* v = S.top();
		S.pop();
		if(visited.find(v) != visited.end()) {
			visited.insert(v);
			if(v->isLeaf())
				n++;
			for(vector<PT>::const_iterator it = v->children.begin(); it != v->children.end(); ++it) {
				const PT* p = &*it;
				S.push(p);
			}
		}
	}
	return n;
}

int PhyloTree::readNiwickTree(const string& treefn) {
	namespace qi = boost::spirit::qi;

	ifstream in(treefn.c_str());
	if(!in.is_open()) {
		cerr << "Error: Could not open input file: " << treefn << endl;
		return -1;
	}

	EGriceLab::newick_grammar<std::string::const_iterator> grammar;
	EGriceLab::PhyloTree tree;

	string content;

	/* copy the whole content in input */
	in.unsetf(std::ios::skipws);
	std::copy(std::istream_iterator<char>(in), std::istream_iterator<char>(), std::back_inserter(content));

	string::const_iterator iter = content.begin();
	string::const_iterator end = content.end();
	// parse
	bool result = qi::phrase_parse(iter, end, grammar, qi::space, tree);

	if(result && iter == end)
		return numNodes();
	else {
		cout << "Parsing Newick tree failed." << endl;
		return -1;
	}
}

ostream& PhyloTree::print(ostream& out) const {
	bool first = true;
	if(!children.empty()) {
		out << '(';
		for(std::vector<PT>::const_iterator it = children.begin(); it != children.end(); ++it) {
			out << (first ? "" : ",");
			it->print(out);
			first = false;
		}
		out << ')';
	}
	if(!name.empty())
		out << name;
	if(length > 0)
		out << ':' << length;

	return out;
}

int PhyloTree::readSeqFromMSA(const MSA* msa) {
	/* check uniqueness of seq-names */
	map<string, unsigned> nameIdx;
	for(unsigned i = 0; i != msa->getNumSeq(); ++i) {
		string name = msa->seqNameAt(i);
		if(nameIdx.find(name) != nameIdx.end()) {
			cerr << "Non-unique seq name " << name << " found in your MSA data " << msa->getName() << endl;
			return -1;
		}
		else
			nameIdx[name] = i;
	}

}

} /* namespace EGriceLab */

