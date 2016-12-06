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
		bool isRooted() const {
			return parent != NULL;
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
	PhyloTreeUnrooted() : L(0) {  }

	/** Construct a PTUnrooted from a Newick Tree */
	PhyloTreeUnrooted(const NewickTree& ntree);

	/* member methods */

	/** Get the number of nodes of this tree using Dfs search */
	long numNodes() const {
		return id2node.size();
	}

	/** get number of aligned sites */
	int getNumAlignSites() const {
		return L;
	}

	/** Get root node */
	const PTUNodePtr& getRoot() const {
		return root;
	}

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
	PTUNodePtr setRoot(PTUNodePtr newRoot);

	/**
	 * test whether the cost (message) of node u -> v has been evaluated
	 */
	bool isEvaluated(const PTUNodePtr u, const PTUNodePtr v) const {
		CostMap::const_iterator result = node2cost.find(u);
		return result != node2cost.end() &&
				result->second.find(v) != result->second.end();
	}

	/**
	 * reset the evaluated cost
	 * @return  actual number of cost removed
	 */
	CostMap::size_type resetCost(const PTUNodePtr u, const PTUNodePtr v) {
		CostMap::iterator result = node2cost.find(u);
		if(result == node2cost.end()) /* u not exists */
			return 0;
		return result->second.erase(v);
	}

	/**
	 * evaluate the likelihood at node node by treating its the root of this PTUnrooted
	 * @param node  new root of this tree
	 * @param model  DNA substitution model, must be a time-reversible model
	 * @return  overall cost of observing this tree at this root
	 */
	double evaluate(const PTUNodePtr node, const DNASubModel& model) const;

private:
	int L; /* number of aligned sites */

	PTUNodePtr root; /* root node of this tree */
	vector<PTUNodePtr> id2node; /* indexed tree nodes */

	BranchLenMap id2length; /* branch length index storing all L*L pairs */
	CostMap node2cost; /* cost index between two neighboring nodes */


public:
	/* static fields */
	static const double MIN_EXPONENT;
};

} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREEUNROOTED_H_ */
