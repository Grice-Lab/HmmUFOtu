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
#include <vector>
#include <map>
#include <stdexcept>
#include <cstddef>
#include <eigen3/Eigen/Dense>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include "StringUtils.h"
#include "SeqCommons.h"
#include "MSA.h"

namespace EGriceLab {

using std::string;
using std::set;
using Eigen::Matrix4Xd;

class PhyloTree {
	/* nested types and enums */
private:
	enum NodeType { P /* parent */, L /* left */, R /* right */ };
	class PhyloTreeNode;

	typedef PhyloTreeNode PTNode;

	class PhyloTreeNode {
	public:
		/* constructors */
		/* Default constructor */
		PhyloTreeNode() : length(0) { }

		/* construct a node with a given DigitalSeq */
		explicit PhyloTreeNode(const DigitalSeq& seq) : seq(seq), length(0) { }

		/* construct a node with a given PrimarySeq */
		explicit PhyloTreeNode(const PrimarySeq& seq) : seq(seq), length(0) { }

		virtual ~PhyloTreeNode() { }

		/* Member methods */
		/** test parent of this node */
		bool hasParent() const {
			return neighbors.find(P) != neighbors.end();
		}

		/** test Lchild of this node */
		bool hasLChild() const {
			return neighbors.find(L) != neighbors.end();
		}

		/** test Rchild of this node */
		bool hasRChild() const {
			return neighbors.find(R) != neighbors.end();
		}

		/** get parent of this node */
		const PTNode& getParent() const {
			return neighbors.find(P)->second;
		}

		/** get Lchild of this node */
		const PTNode& getLChild() const {
			return neighbors.find(L)->second;
		}

		/** get Rchild of this node */
		const PTNode& getRChild() const {
			return neighbors.find(R)->second;
		}

		/** test whether this node is a root node */
		bool isRoot() const;

		/** test whether this node is a leaf node */
		bool isLeaf() const;

		/** test whether this node is a tip node */
		bool isTip() const;

		DigitalSeq seq;
		map<NodeType, PTNode> neighbors;
		double length; /* branch length of this node (to its parent) */

		Matrix4Xd cost; /* cost (negative log liklihood) of observing this sequence given the model and the tree */
	};

	/* constructors */
public:
	/** Default constructor, only assign the alphabet for the entire tree */
	PhyloTree() : abc(SeqCommons::nuclAbc) { }

	/* member methods */
	/** test whether this tree is properly rooted */
	bool isRooted() const;

	/** test whether this tree is leaf rooted */
	bool isLeafRooted() const;

	/** test whether this tree is (arbitrarily) internally rooted */
	bool isInternallyRooted() const;

	/** Get the number of aligned sites of this tree */
	int alnSites() const;

	/** Get the number of nodes of this tree */
	int numNodes() const;

	/** test whether two nodes are siblings */
	bool isSibling(const PhyloTreeNode& node1, const PhyloTreeNode& node2) const;

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

	/**
	 * Read a tree file and MSA into this object
	 * @param treefn  tree filename
	 * @param format  tree file format
	 * @param msa  Multiple Sequence Alignment of this tree
	 * @throw illegal_argument exception if is not a supported file format
	 */
	int readTree(const string& treefn, const string& format, const MSA* msa);

	int readTreeNewick(const string& treefn, const MSA* msa);

private:
	void clear(); /* clear all nodes associated with this PhyloTree */

	PhyloTreeNode root;
	const DegenAlphabet* abc;
	//const DNASubModel* dnaModel;
	//vector<PhyloTreeNode*> nodes;
};

inline bool PhyloTree::PhyloTreeNode::isRoot() const {
	return neighbors.size() == 3 /* an internally rooted tree */
			|| neighbors.size() == 1 /* a leaf root */;
}

inline bool PhyloTree::PhyloTreeNode::isLeaf() const {
	return hasLChild() && hasRChild();
}

inline bool PhyloTree::PhyloTreeNode::isTip() const {
	return !isLeaf() && getLChild().isLeaf() && getRChild().isLeaf();
}

inline bool PhyloTree::isRooted() const {
	return root.isRoot();
}

inline bool PhyloTree::isLeafRooted() const {
	return root.neighbors.size() == 1;
}

inline bool PhyloTree::isInternallyRooted() const {
	return root.neighbors.size() == 3;
}

inline int PhyloTree::alnSites() const {
	return root.seq.length();
}

inline int PhyloTree::numNodes() const {
	return dfsNodes().size();
}

inline PhyloTree::~PhyloTree() {
	clear();
}

inline bool PhyloTree::isSibling(const PhyloTreeNode& node1, const PhyloTreeNode& node2) const {
	assert(isRooted()); // only rooted tree can test siblings
	return !node1.isRoot() && !node2.isRoot() && &node1.getParent() == &node2.getParent();
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

inline int PhyloTree::readTree(const string& treefn, const string& format, const MSA* msa) {
	if(StringUtils::toLower(format) == "newick")
		return readTreeNewick(treefn, msa);
	else
		throw invalid_argument("Unsupported tree file format '" + format + "'");
}

} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREE_H_ */
