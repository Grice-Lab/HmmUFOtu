/*
 * PhyloTreeNode.h
 *	An binary Phylogenic tree class
 *	the tree is rooted if the root node has a null parent, or unrooted if the parent is actually the 3'rd child
 *  Created on: Mar 25, 2016
 *      Author: zhengqi
 */

#ifndef SRC_PHYLOTREE_H_
#define SRC_PHYLOTREE_H_

#include <string>
#include <set>
#include <cstddef>
#include <eigen3/Eigen/Dense>
#include "SeqCommons.h"

namespace EGriceLab {

using std::string;
using std::set;
using Eigen::Matrix4Xd;

class PhyloTree {
	/* nested types and enums */
private:
	class PhyloTreeNode {
	public:
		/* constructors */
		/* Default constructor */
		PhyloTreeNode() : parent(NULL), childL(NULL), childR(NULL),
		parentDist(0), childLDist(0), childRDist(0) { }

		/* construct a node with a given DigitalSeq */
		explicit PhyloTreeNode(const DigitalSeq& seq) :
				seq(seq), parent(NULL), childL(NULL), childR(NULL),
				parentDist(0), childLDist(0), childRDist(0) { }

		/* construct a node with a given PrimarySeq */
		explicit PhyloTreeNode(const PrimarySeq& seq) :
				seq(seq), parent(NULL), childL(NULL), childR(NULL),
				parentDist(0), childLDist(0), childRDist(0) { }

		virtual ~PhyloTreeNode() { }

		/* Member methods */
		/** test whether this node is a root node */
		bool isRoot() const;

		/** test whether this node is a leaf node */
		bool isLeaf() const;

		/** test whether this node is a tip node */
		bool isTip() const;

		DigitalSeq seq;

		PhyloTreeNode* parent;
		PhyloTreeNode* childL;
		PhyloTreeNode* childR;

		double parentDist;
		double childLDist;
		double childRDist;

		Matrix4Xd logLik; /* logLiklihood of observing this sequence given the modeland the tree */
		Matrix4Xd childLPr; /* stored probability matrix for childL */
		Matrix4Xd childRPr; /* stored probability matrix for childR */
		Matrix4Xd parentPr; /* stored probability matrix for parent, only for unrooted tree 'root' */
	};

public:
	/* constructors */
	/** Default constructor*/
	PhyloTree() : root(NULL), abc(SeqCommons::nuclAbc) { }

	/** big-3: Copy constructor */
	PhyloTree(const PhyloTree& tOther);

	/** big-3: copy assignment operator */
	PhyloTree& operator=(PhyloTree tOther); // pass by value intended

	/** big-3: destructor */
	virtual ~PhyloTree();

	/* member methods */
	/** test whether this tree is (arbitrary) rooted */
	bool isRooted() const;
	/** Get the number of aligned sites of this tree */
	int numSites() const;

	/** Get the number of nodes of this tree */
	int numNodes() const;

	/** test whether two nodes are siblings */
	bool isSibling(const PhyloTreeNode* node1, const PhyloTreeNode* node2) const;

	/** Get the DFS set of nodes of a sub tree */
	set<const PhyloTree::PhyloTreeNode*> dfsNodes(const PhyloTree::PhyloTreeNode* node) const;

	/** Get the DFS set of this tree */
	set<const PhyloTree::PhyloTreeNode*> dfsNodes() const;

	/** Get the DFS set of leaf nodes of a sub tree */
	set<const PhyloTree::PhyloTreeNode*> dfsLeaves(const PhyloTree::PhyloTreeNode* node) const;

	/** Get the DFS set of leaf nodes of this tree */
	set<const PhyloTree::PhyloTreeNode*> dfsLeaves() const;

	/** Get the DFS set of tip nodes of a sub tree */
	set<const PhyloTree::PhyloTreeNode*> dfsTips(const PhyloTree::PhyloTreeNode* node) const;

	/** Get the DFS set of tip nodes of this tree */
	set<const PhyloTree::PhyloTreeNode*> dfsTips() const;

	/** Get the parent of a given node */
	set<const PhyloTree::PhyloTreeNode*> parent(const PhyloTree::PhyloTreeNode* node) const {
		return node->parent;
	}

	/** Get the children set of a given node */
	set<const PhyloTree::PhyloTreeNode*> children(const PhyloTree::PhyloTreeNode* node) const;

	/** Get the ancestor set of a given node */
	set<const PhyloTree::PhyloTreeNode*> ancestors(const PhyloTree::PhyloTreeNode* node) const;

	/** Get the offspring set of a given node */
	set<const PhyloTree::PhyloTreeNode*> offsprings(const PhyloTree::PhyloTreeNode* node) const;

	/** Get the sibling of given node of this tree */
	const PhyloTreeNode* getSibling(const PhyloTreeNode*) const;
	PhyloTreeNode* getSibling(const PhyloTreeNode*);

private:
	/* private methods */
	void swap(PhyloTree& other); // no-throw

private:
	PhyloTreeNode* root;
	const DegenAlphabet* abc;
	//const DNASubModel* dnaModel;
	//vector<PhyloTreeNode*> nodes;
};

inline bool PhyloTree::PhyloTreeNode::isRoot() const {
	return parent == NULL;
}

inline bool PhyloTree::PhyloTreeNode::isLeaf() const {
	return childL == NULL && childR == NULL;
}

inline bool PhyloTree::PhyloTreeNode::isTip() const {
	return !isLeaf() && childL->isLeaf() && childR->isLeaf();
}

inline void PhyloTree::swap(PhyloTree& other) {
	using std::swap;
	swap(root, other.root);
	swap(abc, other.abc);
}

inline bool PhyloTree::isRooted() const {
	return root->isRoot();
}

inline int PhyloTree::numSites() const {
	return root->seq.length();
}

inline int PhyloTree::numNodes() const {
	return dfsNodes().size();
}

inline bool PhyloTree::isSibling(const PhyloTreeNode* node1, const PhyloTreeNode* node2) const {
	assert(isRooted()); // only rooted tree can test siblings
	return !node1->isRoot() && !node2->isRoot() && node1->parent == node2->parent;
}

inline set<const PhyloTree::PhyloTreeNode*> PhyloTree::children(const PhyloTree::PhyloTreeNode* node) const {
	set<const PhyloTree::PhyloTreeNode*> children;
	children.insert(node->childL);
	children.insert(node->childR);
	return children;
}

inline set<const PhyloTree::PhyloTreeNode*> PhyloTree::ancestors(const PhyloTree::PhyloTreeNode* node) const {
	set<const PhyloTree::PhyloTreeNode*> ancestors;
	while(!node->isRoot()) {
		ancestors.insert(node->parent);
		node = node->parent;
	}
	return ancestors;
}

inline set<const PhyloTree::PhyloTreeNode*> PhyloTree::offsprings(const PhyloTree::PhyloTreeNode* node) const {
	set<const PhyloTree::PhyloTreeNode*>& offsprings = dfsNodes(node);
	offsprings.erase(node); /* remove itself */
	return offsprings;
}

inline const PhyloTree::PhyloTreeNode* PhyloTree::getSibling(const PhyloTreeNode* node) const {
	assert(isRooted()); /* only rooted tree can find siblings */
	if(node->isRoot())
		return NULL;
	return node->parent->childL == node ? node->parent->childR : node->parent->childL;
}

inline PhyloTree::PhyloTreeNode* PhyloTree::getSibling(const PhyloTreeNode* node) {
	assert(isRooted()); /* only rooted tree can find siblings */
	if(node->isRoot())
		return NULL;
	return node->parent->childL == node ? node->parent->childR : node->parent->childL;
}

} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREE_H_ */
