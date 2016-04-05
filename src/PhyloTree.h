/*
 * PhyloTreeNode.h
 *
 *  Created on: Mar 25, 2016
 *      Author: zhengqi
 */

#ifndef SRC_PHYLOTREE_H_
#define SRC_PHYLOTREE_H_

#include <string>
#include <vector>
#include <cstddef>

namespace EGriceLab {

using std::string;
using std::vector;

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
	};

public:
	/* constructors */
	/** Default constructor*/
	PhyloTree() : root(NULL) { }

	/** big-3: Copy constructor */
	PhyloTree(const PhyloTree& tOther);

	/** big-3: copy assignment operator */
	PhyloTree& operator=(PhyloTree tOther); // pass by value intended

	/** big-3: destructor */
	virtual ~PhyloTree();

private:
	/* private methods */
	void swap(PhyloTree& other); // no-throw

private:
	PhyloTreeNode* root;
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
}

} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREE_H_ */
