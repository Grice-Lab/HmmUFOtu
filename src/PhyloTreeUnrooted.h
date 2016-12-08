/*
 * PhyloTreeUnrooted.h
 *  An Unrooted Phylogenic Tree (PTUnrooted)
 *  A PTUnrooted can be evaluated from any node as its root and yields same cost
 *  as long as using a time-reversible DNA substitution model
 *  Internal Tree nodes are number indexed from 0 to N-1
 *  Created on: Dec 1, 2016
 *      Author: zhengqi
 */

#ifndef SRC_PHYLOTREEUNROOTED_H_
#define SRC_PHYLOTREEUNROOTED_H_

#include <string>
#include <vector>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <cstddef>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include "HmmUFOtuConst.h"
#include "StringUtils.h"
#include "SeqCommons.h"
#include "DigitalSeq.h"
#include "NewickTree.h"
#include "MSA.h"
#include "DNASubModel.h"

namespace EGriceLab {
using std::string;
using std::vector;
using std::map;
using std::istream;
using std::ostream;
using Eigen::Matrix4Xd;
using Eigen::RowVectorXd;
using boost::shared_ptr;
using boost::unordered_map;

class PhyloTreeUnrooted; /* forward declaration */

typedef PhyloTreeUnrooted PTUnrooted;

class PhyloTreeUnrooted {
public:
	/* nested types and enums */

	struct PhyloTreeUnrootedNode;
	typedef PTUnrooted::PhyloTreeUnrootedNode PTUNode;

	typedef shared_ptr<PTUNode> PTUNodePtr; /* use boost shared_ptr to hold node pointers */


	/**
	 * A PTUnrooed node that stores its basic information and neighbors
	 */
	struct PhyloTreeUnrootedNode {
		/* constructors */
		/**
		 * Default constructor, do nothing
		 */
		PhyloTreeUnrootedNode() : id(0), annoDist(0) {	}

		/**
		 * Construct a PTUNode with a given name and id
		 */
		explicit PhyloTreeUnrootedNode(long id, const string& name = "")
		: id(id), name(name), annoDist(0) {  }

		/**
		 * Construct a PTUNode with a given id, name and sequence
		 */
		PhyloTreeUnrootedNode(long id, const string& name, const DigitalSeq& seq)
		: id(id), name(name), seq(seq), annoDist(0)
		{ }

		/* Member methods */

		/** test whether this node is named */
		bool isNamed() const {
			return !name.empty();
		}

		/** test whether this is a leave node */
		bool isLeaf() const {
			return neighbors.size() == 1;
		}

		/** test whether this is an internal node */
		bool isInternal() const {
			return neighbors.size() > 1;
		}

		/** test whether this subtree is rooted */
		bool isRoot() const {
			return parent == NULL;
		}

		/* member fields */
		long id; /* a unique id for each node */
		string name; /* node name, need to be unique for database loading */
		DigitalSeq seq; /* sequence of this node */
		vector<PTUNodePtr> neighbors; /* pointers to neighbors */
		PTUNodePtr parent; /* pointer to parent node, set to null on default */

		string anno;
		double annoDist;
	};

	typedef unordered_map<PTUNodePtr, unordered_map<PTUNodePtr, Matrix4Xd> > CostMap;
	typedef unordered_map<PTUNodePtr, unordered_map<PTUNodePtr, double> > BranchLenMap;

	/* constructors */
	/** Default constructor, do nothing */
	PhyloTreeUnrooted() : csLen(0) {  }

	/** Construct a PTUnrooted from a Newick Tree */
	PhyloTreeUnrooted(const NewickTree& ntree);

	/* member methods */

	/** Get the number of nodes of this tree */
	size_t numNodes() const {
		return id2node.size();
	}

	/** Get number of leaves in this tree */
	size_t numLeaves() const;

	/** get number of aligned sites */
	int getNumAlignSites() const {
		return csLen;
	}

	/** Get root node */
	const PTUNodePtr& getRoot() const {
		return root;
	}

	/* get all nodes */
	std::vector<PTUNodePtr> getNodes() const {
		return id2node;
	}

	/* get node i */
	PTUNodePtr getNode(std::vector<PTUNodePtr>::size_type i) const {
		return id2node[i];
	}

	/**
	 * get branch length from u -> v
	 * @return  -1 if not exists
	 */
	double getBranchLength(const PTUNodePtr& u, const PTUNodePtr& v) const;

	/** Load sequences from MSA into this this */
	long loadMSA(const MSA& msa);

	/**
	 * save this object to an output in binary format
	 */
	ostream& save(ostream& out) const;

	/** load a binary file into this object */
	istream& load(string& in);

	/**
	 * set tree root at given node, return the old node
	 */
	PTUNodePtr setRoot(const PTUNodePtr& newRoot);

	/**
	 * test whether the cost (message) of node u -> v of site j has been evaluated
	 */
	bool isEvaluated(const PTUNodePtr u, const PTUNodePtr v, int j) const {
		return (node2cost.find(u)->second.find(v)->second.col(j).array() != inf).any();
	}

	/**
	 * initiate the cached incoming cost of between every node u and every neighbor v
	 */
	void initInCost();

	/**
	 * initiate the leaf cost given a DNA Sub-model for faster computation
	 * @param model  any Time-reversible DNA substitution model
	 */
	void initLeafCost(const DNASubModel& model);

	/**
	 * reset the cached cost of edge u->v
	 */
	void resetCost(const PTUNodePtr u, const PTUNodePtr v) {
		node2cost[u][v].setConstant(inf);
	}

	/**
	 * reset the cached cost of every node
	 */
	void resetCost();

	/**
	 * evaluate the likelihood at given node by treating it as the root of this PTUnrooted
	 * @param node  new root of this tree
	 * @param model  DNA substitution model, must be a time-reversible model
	 * @return  cost matrix of observing this tree at this root
	 */
	Matrix4Xd evaluate(const PTUNodePtr& node, const DNASubModel& model);

	/**
	 * evaluate the likelihood of a given node and site by treating it as the root of this PTUnrooted
	 * @param node  new root of this tree
	 * @param j  the jth aligned site
	 * @param model  DNA substitution model, must be a time-reversible model
	 * @return  cost vector of observing this tree at the jth site at this root
	 */
	Vector4d evaluate(const PTUNodePtr& node, int j, const DNASubModel& model);

	/**
	 * calculate the entire tree cost given a DNA sub model
	 * evaluate the tree if necessary
	 */
	double treeCost(const DNASubModel& model);

	/**
	 * calculate the entire tree cost at given region using given a DNA sub model
	 * both start and end are 0-based inclusive
	 */
	double treeCost(int start, int end, const DNASubModel& model);

	/**
	 * calculate the entire tree cost at j-th site given a DNA sub model
	 * evaluate the tree if neccessary
	 */
	double treeCost(int j, const DNASubModel& model);

	/**
	 * write this PTUnrooted tree structure into output in given format
	 */
	ostream& writeTree(ostream& out, string format = "newick") const;

	/**
	 * write this PTUnrooted tree structure with given root recursively into output in newick format
	 */
	ostream& writeTreeNewick(ostream& out, const PTUNodePtr& node) const;

	/* static methods */
	/**
	 * test whether p is parent of c
	 */
	static bool isParent(const PTUNodePtr& p, const PTUNodePtr& c) {
		return c != NULL && c->parent == p;
	}

	static bool isChild(const PTUNodePtr& c, const PTUNodePtr& p) {
		return isParent(p, c);
	}

	/* member fields */
private:
	int csLen; /* number of aligned sites */

	PTUNodePtr root; /* root node of this tree */
	vector<PTUNodePtr> id2node; /* indexed tree nodes */

	BranchLenMap node2length; /* branch length index storing all L*L pairs */
	CostMap node2cost; /* cached cost message sending from u -> v, before conjugating into the Pr(v) of the branch-length */
	Matrix4Xd leafCost; /* cached 4 X 5 leaf cost matrix,
						with each column the pre-computed cost of observing A, C, G, T or - at any given site */

public:
	/* static fields */
	static const double MAX_COST_EXP;
};

inline size_t PTUnrooted::numLeaves() const {
	size_t nLeaves = 0;
	for(vector<PTUNodePtr>::const_iterator nodeIt = id2node.begin(); nodeIt != id2node.end(); ++nodeIt)
		if((*nodeIt)->isLeaf())
			nLeaves++;
	return nLeaves;
}

inline Matrix4Xd PTUnrooted::evaluate(const PTUNodePtr& node, const DNASubModel& model) {
	Matrix4Xd costMat(4, csLen);
	for(int j = 0; j < csLen; ++j) {
		std::cerr << "Evaluating at site " << j << std::endl;
		costMat.col(j) = evaluate(node, j, model);
	}
	return costMat;
}

inline ostream& PTUnrooted::writeTree(ostream& out, string format) const {
	StringUtils::toLower(format);
	if(format == "newick")
		return writeTreeNewick(out, root) << ";";
	else {
		std::cerr << "Cannot write PTUnrooted tree, unknown tree format " << format << std::endl;
		out.setstate(std::ios_base::failbit);
		return out;
	}
}

inline double PTUnrooted::getBranchLength(const PTUNodePtr& u, const PTUNodePtr& v) const {
	BranchLenMap::const_iterator resultOuter = node2length.find(u);
	if(resultOuter == node2length.end())
		return -1;
	else {
		boost::unordered_map<PTUNodePtr, double>::const_iterator resultInner = resultOuter->second.find(v);
		return resultInner == resultOuter->second.end() ? -1 : resultInner->second;
	}
}

inline double PTUnrooted::treeCost(int start, int end, const DNASubModel& model) {
	double cost = 0;
	for(int j = start; j <= end; ++j)
		cost += treeCost(j, model);
	return cost;
}

inline double PTUnrooted::treeCost(const DNASubModel& model) {
	return treeCost(0, csLen - 1, model);
}


} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREEUNROOTED_H_ */
