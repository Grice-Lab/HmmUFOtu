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

PhyloTree::PhyloTree(const PhyloTree& tOther) : root(NULL), abc(tOther.abc) {
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

set<const PhyloTree::PhyloTreeNode*> PhyloTree::dfsNodes(const PhyloTree::PhyloTreeNode* node) const {
	assert(isRooted());
	/* Do a DFS explore of the tree */
	stack<const PhyloTreeNode*> S; // using a stack to do non-recursive DFS
	set<const PhyloTreeNode*> visited; // using a set to determine whether this node has been visited
	S.push(node);
	while(!S.empty()) {
		PhyloTreeNode* v = S.pop();
		if(!visited.find(v)) { // this node has not been visited
			visited.insert(v);
			if(!v->isLeaf()) {
				S.push(v->childL);
				S.push(v->childR);
			}
		}
	}
	return visited;
}

set<const PhyloTree::PhyloTreeNode*> PhyloTree::dfsNodes() const {
	return dfsNodes(root);
}

set<const PhyloTree::PhyloTreeNode*> PhyloTree::dfsLeaves(const PhyloTree::PhyloTreeNode* node) const {
	/* Do a DFS explore of the tree */
	stack<const PhyloTreeNode*> S; // using a stack to do non-recursive DFS
	set<const PhyloTreeNode*> visited; // using a set to determine whether this node has been visited
	set<const PhyloTreeNode*> leaves;
	S.push(node);
	while(!S.empty()) {
		PhyloTreeNode* v = S.pop();
		if(!visited.find(v)) { // this node has not been visited
			visited.insert(v);
			if(!v->isLeaf()) {
				S.push(v->childL);
				S.push(v->childR);
			}
			else
				leaves.insert(v);
		}
	}
	return leaves;
}

set<const PhyloTree::PhyloTreeNode*> PhyloTree::dfsLeaves() const {
	return dfsLeaves(root);
}

set<const PhyloTree::PhyloTreeNode*> PhyloTree::dfsTips(const PhyloTree::PhyloTreeNode* node) const {
	/* Do a DFS explore of the tree */
	stack<const PhyloTreeNode*> S; // using a stack to do non-recursive DFS
	set<const PhyloTreeNode*> visited; // using a set to determine whether this node has been visited
	set<const PhyloTreeNode*> tips;
	S.push(node);
	while(!S.empty()) {
		PhyloTreeNode* v = S.pop();
		if(!visited.find(v)) { // this node has not been visited
			visited.insert(v);
			if(!v->isLeaf()) {
				S.push(v->childL);
				S.push(v->childR);
			}
			if(v->isTip())
				tips.insert(v);
		}
	}
	return tips;
}

set<const PhyloTree::PhyloTreeNode*> PhyloTree::dfsTips() const {
	return dfsTips(root);
}

} /* namespace EGriceLab */

