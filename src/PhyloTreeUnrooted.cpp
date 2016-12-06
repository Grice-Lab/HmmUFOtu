/*
 * PhyloTreeUnrooted.cpp
 *
 *  Created on: Dec 1, 2016
 *      Author: zhengqi
 */

#include <stack>
#include <boost/unordered_set.hpp>
#include "PhyloTreeUnrooted.h"

namespace EGriceLab {
using namespace std;
using namespace EGriceLab;

const double PhyloTreeUnrooted::MIN_EXPONENT = std::numeric_limits<float>::min_exponent;

PhyloTreeUnrooted::PhyloTreeUnrooted(const NewickTree& ntree) : L(0) {
	/* DFS of the NewickTree */
	boost::unordered_set<const NT*> visited;
	stack<const NT*> S;
	long id = 0; /* id start from 0 */

	S.push(&ntree);
	while(!S.empty()) {
		const NT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* construct this PTUNode */
			PTUNodePtr u(new PTUNode(id++, v->name));
			id2node.push_back(u);

			/* set root, if this is the first encounted node */
			if(root == NULL)
				root = u;

			/* check each child of u(v) */
			for(vector<NT>::const_iterator childIt = v->children.begin(); childIt != v->children.end(); ++childIt) {
				/* construct a PTUNode based on this child */
				PTUNodePtr child = PTUNodePtr(new PTUNode(id++, childIt->name));
				/* update the child */
				child->parent = u;
				child->neighbors.push_back(u);
				/* update the parent */
				u->neighbors.push_back(child);

				const NT* p = &*childIt;
				S.push(p);
			}
		}
	}
}

long PhyloTreeUnrooted::loadMSA(const MSA& msa) {
	const DegenAlphabet* abc = msa.getAbc();
	if(abc != SeqCommons::nuclAbc) {
		cerr << "PhyloTreeUnrooted can only read in MSA in dna alphabet" << endl;
		return -1;
	}
	const unsigned numSeq = msa.getNumSeq();
	const unsigned csLen = msa.getCSLen();
	assert(numNodes() >= numSeq);

	/* check uniqueness of seq names in msa */
	map<string, unsigned> nameIdx;
	for(unsigned i = 0; i != numSeq; ++i) {
		string name = msa.seqNameAt(i);
		if(nameIdx.find(name) != nameIdx.end()) {
			cerr << "Non-unique seq name " << name << " found in your MSA data " << msa.getName() << endl;
			return -1;
		}
		else {
			nameIdx[name] = i;
		}
	}
	/* assign seq to each nodes of the tree, exit early if cannot find */
	for(vector<PTUNodePtr>::size_type id = 0; id != id2node.size(); ++id) {
		const PTUNodePtr& node = id2node[id];
		assert(id == node->id);

		map<string, unsigned>::const_iterator result = nameIdx.find(node->name);
		if(result == nameIdx.end()) { /* this name cannot be found in the msa */
			cerr << "Tree node " << node->name << " does not exist in the MSA data" << endl;
			return id; /* id can be 0 */
		}
		node->seq = msa.dsAt(nameIdx[node->name]);
	}
	return id2node.size();
}

PhyloTreeUnrooted::PTUNodePtr PhyloTreeUnrooted::setRoot(PTUNodePtr newRoot) {
	if(newRoot == NULL || newRoot == root) /* no need to set */
		return root;

	newRoot->parent = NULL; // root has no parent
	/* DFS of this tree starting from newRoot */
	boost::unordered_set<PTUNodePtr> visited;
	stack<PTUNodePtr> S;

	S.push(newRoot);
	while(!S.empty()) {
		PTUNodePtr v = S.top();
		S.pop();

		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);

			/* check each child of v */
			for(vector<PTUNodePtr>::iterator neighborIt = v->neighbors.begin(); neighborIt != v->neighbors.end(); ++neighborIt) {
				if(*neighborIt == v) // this is neighbor is its parent, ignore
					continue;
				/* update the child's parent */
				(*neighborIt)->parent = v;
				S.push(*neighborIt);
			}
		}
	}
	PTUNodePtr oldRoot = root;
	root = newRoot;
	return oldRoot;
}

} /* namespace EGriceLab */
