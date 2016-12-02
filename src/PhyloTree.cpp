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
#include <limits>
#include "PhyloTree.h"

namespace EGriceLab {
using namespace std;

const double PhyloTree::MIN_EXPONENT = std::numeric_limits<float>::min_exponent;

int PhyloTree::numNodes() const {
	set<const PT*> visited;
	stack<const PT*> S;

	S.push(this);
	while(!S.empty()) {
		const PT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
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
		if(visited.find(v) == visited.end()) { /* not visited before */
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

set<const PT*> PhyloTree::leafNodes() const {
	set<const PT*> visited;
	stack<const PT*> S;
	set<const PT*> leaves;

	S.push(this);
	while(!S.empty()) {
		const PT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			if(v->isLeaf())
				leaves.insert(v);
			for(vector<PT>::const_iterator it = v->children.begin(); it != v->children.end(); ++it) {
				const PT* p = &*it;
				S.push(p);
			}
		}
	}
	return leaves;
}

const PT* PhyloTree::firstLeaf() const {
	const PT* node = this;
	while(!node->isLeaf())
		node = &(node->children.front());
	return node;
}


int PhyloTree::readNiwickTree(const string& treefn) {
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
	if(!name.empty()) {
		if(name.find('\'') == string::npos) // name doesn't contains quotes
			out << name;
		else
			out << "'" << name << "'";
	}
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
		else {
			nameIdx[name] = i;
		}
	}
	/* ass Digital seq to each nodes of the tree using DFS */
	int n = 0;
	set<PT*> visited;
	stack<PT*> S;

	S.push(this);
	while(!S.empty()) {
		PT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* process this node */
			v->cost.resize(4, msa->getCSLen()); /* initiate the cost */
			v->cost.setConstant(inf);
			/* initiate the scale */
			v->scale.resize(msa->getCSLen());
			v->scale.setZero();
			if(nameIdx.find(v->name) != nameIdx.end()) { // this node exists in MSA
				v->seq = msa->dsAt(nameIdx[v->name]);
				n++;
			}
			else /* assign an all-gap seq for this node */
				v->seq = DigitalSeq(msa->getAbc(), v->name, string(msa->getCSLen(), msa->getAbc()->getGap()[0]));

			for(vector<PT>::iterator it = v->children.begin(); it != v->children.end(); ++it) {
				PT* p = &*it;
				S.push(p);
			}
		}
	}
	return n;
}

long PhyloTree::setIDandParent() {
	long currID = 0;
	/* browse the tree in Depth-first order */
	set<PT*> visited;
	stack<PT*> S;

	S.push(this);
	while(!S.empty()) {
		PT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* set the id */
			v->id = currID++;
			/* set the parent of its children */
			for(vector<PT>::iterator it = v->children.begin(); it != v->children.end(); ++it) {
				it->parent = v;
				PT* p = &*it;
				S.push(p);
			}
		}
	}
	return currID;
}

long PhyloTree::updateParent() {
	/* browse the tree in Depth-first order */
	set<PT*> visited;
	stack<PT*> S;

	S.push(this);
	while(!S.empty()) {
		PT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* set the parent of its children */
			for(vector<PT>::iterator it = v->children.begin(); it != v->children.end(); ++it) {
				it->parent = v;
				PT* p = &*it;
				S.push(p);
			}
		}
	}
}

void PhyloTree::annotate() {
	/* browse the tree in Depth-first order */
	set<PT*> visited;
	stack<PT*> S;

	S.push(this);
	while(!S.empty()) {
		PT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* trace back ancestors to annotate this node */
			for(const PT* p = v; p->parent != NULL; p = p->parent) {
				if(p->isNamed()) {
					v->anno = p->name;
					break;
				}
				v->annoDist += p->length; /* add up length */
			}
			if(v->anno.empty()) /* if cannot find any named ancestor */
				v->anno = "root";

			for(vector<PT>::iterator it = v->children.begin(); it != v->children.end(); ++it) {
				PT* p = &*it;
				S.push(p);
			}
		}
	}
}

} /* namespace EGriceLab */

