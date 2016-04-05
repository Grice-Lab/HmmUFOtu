/*
 * PhyloTreeNode.cpp
 *
 *  Created on: Mar 25, 2016
 *      Author: zhengqi
 */

#include "PhyloTree.h"
#include <stack>
#include <set>

namespace EGriceLab {
using namespace std;

PhyloTree::PhyloTree(const PhyloTree& tOther) : root(NULL) {
	/* do a DFS copy of the tree edges starting from root */
	stack<PhyloTreeNode*> S; // using a stack to do non-recursive DFS
	set<PhyloTreeNode*> copied; // using a set to determine whether this node has been copied

	S.push(tOther.root);

	while(!S.empty()) {
		PhyloTreeNode* v = S.pop();
		if(!copied.find(v)) { // this node has not been copied
			copied.insert(v);
			/* Make a new copy of this root */
			PhyloTreeNode* newV = new PhyloTreeNode(*v);
			// reset its parent and children
			newV->parent = newV->childL = newV->childR = NULL;
			/* copy childL */
			if(v->childL != NULL) {
				PhyloTreeNode* newChildL = new PhyloTreeNode(v->childL);
				newV->childL = newChildL;
				newChildL->parent = newV;
				newChildL->childL = NULL;
				newChildL->childR = NULL;
				S.push(v->childL);
			}
			/* copy childR */
			if(v->childR != NULL) {
				PhyloTreeNode* newChildR = new PhyloTreeNode(v->childR);
				newV->childR = newChildR;
				newChildR->parent = newV;
				newChildR->childL = NULL;
				newChildR->childR = NULL;
				S.push(v->childR);
			}
		}
	}
}

PhyloTree& PhyloTree::operator =(PhyloTree tOther) { // passed by value
	tOther.swap(*this); // swap the copied temporary parameter
	return *this; // old resources now in *this
}

PhyloTree::~PhyloTree() {
	/* do a DFS deconstruction of the entire tree from root */
	stack<PhyloTreeNode*> S; // using a stack to do non-recursive DFS
	S.push(root);

	while(!S.empty()) {
		PhyloTreeNode* v = S.pop();
		if(v != NULL) { // this node has not been destroyed
			if(v->childL != NULL)
				S.push(v->childL);
			if(v->childR != NULL)
				S.push(v->childR);
			/* destroy this node */
			delete v;
			v = NULL;
		}
	}
}

} /* namespace EGriceLab */

