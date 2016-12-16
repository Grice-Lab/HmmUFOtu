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
#include <cstdlib>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include "HmmUFOtuConst.h"
#include "ProgLog.h"
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
using Eigen::Matrix4d;
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
	typedef shared_ptr<const PTUNode> PTUNodeConstPtr; /* use boost shared_ptr to hold node pointers */

	/**
	 * A PTUnrooed node that stores its basic information and neighbors
	 */
	class PhyloTreeUnrootedNode {
		friend class PhyloTreeUnrooted;

	public:
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
		/* Getters and Setters */
		const string& getAnno() const {
			return anno;
		}

		double getAnnoDist() const {
			return annoDist;
		}

		long getId() const {
			return id;
		}

		const string& getName() const {
			return name;
		}

		const PTUNodePtr& getParent() const {
			return parent;
		}

		const DigitalSeq& getSeq() const {
			return seq;
		}

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

		/** test whether this is a root node */
		bool isRoot() const {
			return parent == NULL;
		}

		/** test whether this node is parent of another node */
		bool isParent(const PTUNodePtr& other) const {
			return other != NULL && this == other->parent.get();
		}

		/** test whether this node is child of another node */
		bool isChild(const PTUNodePtr& other) const {
			return parent == other;
		}

		/**
		 * test whether this is a tip node
		 * all children of a tip must be leaves
		 */
		bool isTip() const {
			if(isLeaf())
				return false;
			for(vector<PTUNodePtr>::const_iterator child = neighbors.begin(); child != neighbors.end(); ++child)
				if(isParent(*child) /* this is really a child */
						&& !(*child)->isLeaf())
					return false;
			return true;
		}

		/**
		 * get children of this node
		 * children nodes are neighbors excluding the parent
		 */
		vector<PTUNodePtr> getChildren() const {
			vector<PTUNodePtr> children (neighbors);
			children.erase(std::remove(children.begin(), children.end(), parent), children.end());
			return children;
		}

		/**
		 * get first child of this node
		 * return NULL if not exists
		 */
		PTUNodePtr firstChild() const {
			for(vector<PTUNodePtr>::const_iterator child = neighbors.begin(); child != neighbors.end(); ++child)
				if(isParent(*child)) // this is really a child
					return *child;
			return NULL;
		}

		/**
		 * get last child of this node
		 * return NULL if not exists
		 */
		PTUNodePtr lastChild() const {
			for(vector<PTUNodePtr>::const_reverse_iterator child = neighbors.rbegin(); child != neighbors.rend(); ++child)
				if(isParent(*child)) // this is really a child
					return *child;
			return NULL;
		}

		/**
		 * get first leaf as an offspring of this node
		 */
		const PTUNode* firstLeaf() const {
			const PTUNode* node = this; /* search from this node */
			while(!node->isLeaf())
				node = node->firstChild().get();
			return node;
		}

		/**
		 * get first leaf as an offspring of this node
		 */
		PTUNode* firstLeaf() {
			PTUNode* node = this; /* search from this node */
			while(!node->isLeaf())
				node = node->firstChild().get();
			return node;
		}

		/**
		 * get the last leaf as an offspring of this node
		 */
		const PTUNode* lastLeaf() const {
			const PTUNode* node = this; /* search from this node */
			while(!node->isLeaf())
				node = node->lastChild().get();
			return node;
		}

		/**
		 * get the last leaf as an offspring of this node
		 */
		PTUNode* lastLeaf() {
			PTUNode* node = this; /* search from this node */
			while(!node->isLeaf())
				node = node->lastChild().get();
			return node;
		}

		/**
		 * get a random leaf as an offspring of this node
		 * you need to call srand() in your main program
		 */
		const PTUNode* randomLeaf() const {
			const PTUNode* node = this; /* search from this node */
			while(!node->isLeaf()) {
				const vector<PTUNodePtr>& children = node->getChildren();
				node = children[rand() % children.size()].get();
			}
			return node;
		}

		/**
		 * get a random leaf as an offspring of this node
		 * you need to call srand() in your main program
		 */
		PTUNode* randomLeaf() {
			PTUNode* node = this; /* search from this node */
			while(!node->isLeaf()) {
				const vector<PTUNodePtr>& children = node->getChildren();
				node = children[rand() % children.size()].get();
			}
			return node;
		}

		int numNeighbors() const {
			return neighbors.size();
		}

	private:
		/**
		 * save this node information but not its edges to a binary file
		 */
		istream& load(istream& in);

		/**
		 * load the node information but not its edges from a binary file
		 */
		ostream& save(ostream& out) const;

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

	/* disable the copy and assignment constructors */
private:
	PhyloTreeUnrooted(const PhyloTreeUnrooted& other);
	PhyloTreeUnrooted& operator=(const PhyloTreeUnrooted& other);

public:

	/* member methods */

	/** Get the number of nodes of this tree */
	size_t numNodes() const {
		return id2node.size();
	}

	size_t numEdges() const;

	/** Get number of leaves in this tree */
	size_t numLeaves() const;

	/** get number of aligned sites */
	int numAlignSites() const {
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

	/**
	 * get branch cost from u->v, before convoluted into the Pr(length)
	 * @return uinitiated matrix if not exists
	 */
	Matrix4Xd getBranchCost(const PTUNodePtr& u, const PTUNodePtr& v) const;

	/** Load sequences from MSA into this this */
	long loadMSA(const MSA& msa);

	/**
	 * save this object to an output in binary format
	 */
	ostream& save(ostream& out) const;

	/** load a binary file into this object */
	istream& load(istream& in);

	/**
	 * set tree root at given node, return the old node
	 */
	PTUNodePtr setRoot(const PTUNodePtr& newRoot);

	/**
	 * set tree root at the ith node, return the old node id
	 */
	size_t setRoot(size_t newRootId) {
		return setRoot(id2node[newRootId])->id;
	}

	/**
	 * test whether the cost (message) of node u -> v of site j has been evaluated
	 */
	bool isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v, int j) const {
		return (node2cost.at(u).at(v).col(j).array() != INVALID_COST).all();
	}

	/**
	 * initiate the cached incoming cost of between every node u and every neighbor v
	 */
	void initInCost();

	/**
	 * initiate the leaf cost to all inf
	 */
	void initLeafCost() {
		leafCost = Matrix4Xd::Constant(4, 5, INVALID_COST);
	}

	/**
	 * initiate the leaf cost given a DNA Sub-model for faster computation
	 * @param model  any Time-reversible DNA substitution model
	 */
	void initLeafCost(const DNASubModel& model);

	/**
	 * reset the cached cost of edge u->v
	 */
	void resetCost(const PTUNodePtr& u, const PTUNodePtr& v) {
		node2cost[u][v].setConstant(INVALID_COST);
	}

	/**
	 * reset the cached cost of every node
	 */
	void resetCost();

	/**
	 * reset the cached leaf cost
	 */
	void resetLeafCost() {
		leafCost.setConstant(INVALID_COST);
	}

	/**
	 * evaluate the likelihood at current root by treating it as the root of this PTUnrooted
	 * @param model  DNA substitution model, must be a time-reversible model
	 */
	void evaluate(const DNASubModel& model) {
		evaluate(root, model);
	}

	/**
	 * evaluate the likelihood at given node by treating it as the root of this PTUnrooted
	 * @param node  new root of this tree
	 * @param model  DNA substitution model, must be a time-reversible model
	 */
	void evaluate(const PTUNodePtr& node, const DNASubModel& model);

	/**
	 * evaluate a given node at the jth site, given a DNA model
	 * @param node  new root of this tree
	 * @param j  the jth aligned site
	 * @param model  DNA substitution model, must be a time-reversible model
	 */
	void evaluate(const PTUNodePtr& node, int j, const DNASubModel& model);

	/**
	 * evaluate the incoming cost u->v at the jth site, given a DNA model
	 * @param node  new root of this tree
	 * @param j  the jth aligned site
	 * @param model  DNA substitution model, must be a time-reversible model
	 * @return  cost vector of observing this tree at the jth site of this edge
	 */
	Vector4d evaluate(const PTUNodePtr& u, const PTUNodePtr& v, int j, const DNASubModel& model);

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

	/**
	 * get a transition dataset for parameter traning of a DNA Substitution model
	 * using one of two well studied method, "Gojobori" or "Goldman"
	 */
	vector<Matrix4d> getModelTransitionSet(string method = "Gojobori") const;

	/**
	 * get training dataset for model parameters training of a DNA Substitution model
	 * using one of two well studied method, "Gojobori" or "Goldman"
	 */
	vector<Matrix4d> getModelTraningSetGojobori() const;

	/**
	 * get training dataset for model parameters training of a DNA Substitution model
	 * using one of two well studied method, "Gojobori" or "Goldman"
	 */
	vector<Matrix4d> getModelTraningSetGoldman() const;

	/**
	 * get estimated base frequency (pi) using this tree
	 */
	Vector4d getModelFreqEst() const;

private:
	/**
	 * load an edge node1->node2 from a binary input
	 */
	istream& loadEdge(istream& in);

	/**
	 * save an edge node1->node2 to a binary output
	 * only the relationship between node IDs are stored
	 */
	ostream& saveEdge(ostream& out, const PTUNodePtr& node1, const PTUNodePtr& node2) const;

	/**
	 * load the edge cost node1->node2 from a binary input
	 */
	istream& loadEdgeCost(istream& in);

	/**
	 * save the edge cost node1->node2 to a binary output
	 */
	ostream& saveEdgeCost(ostream& out, const PTUNodePtr& node, const PTUNodePtr& node2) const;

	/**
	 * load leaf cost from a binary input
	 */
	istream& loadLeafCost(istream& in);

	/**
	 * save leaf cost to a binary outout
	 */
	ostream& saveLeafCost(ostream& out) const;

	/**
	 * load root information from a binary input
	 */
	istream& loadRoot(istream& in);

	/**
	 * save root information to a binary output
	 */
	ostream& saveRoot(ostream& out) const;

public:
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

	/**
	 * test whether a node a tip
	 * all children of a tip must be leaves
	 */
	static bool isTip(const PTUNodePtr& node);

	static PTUNodePtr firstLeaf(PTUNodePtr node);
	static PTUNodePtr lastLeaf(PTUNodePtr node);
	static PTUNodePtr randomLeaf(PTUNodePtr node);

	/* member fields */
private:
	int csLen; /* number of aligned sites */

	PTUNodePtr root; /* root node of this tree */
	vector<PTUNodePtr> id2node; /* indexed tree nodes */

	BranchLenMap node2length; /* branch length index storing edge length */
	CostMap node2cost; /* cached cost message sending from u -> v, before conjugating into the Pr(v) of the branch-length */
	Matrix4Xd leafCost; /* cached 4 X 5 leaf cost matrix,
						with each column the pre-computed cost of observing A, C, G, T or - at any given site */

public:
	/* static fields */
	static const double MAX_COST_EXP;
	static const double INVALID_COST;
};

inline size_t PTUnrooted::numEdges() const {
	size_t nEdges = 0;
	for(vector<PTUNodePtr>::const_iterator nodeIt = id2node.begin(); nodeIt != id2node.end(); ++nodeIt)
		nEdges += (*nodeIt)->numNeighbors();
	return nEdges;
}


inline size_t PTUnrooted::numLeaves() const {
	size_t nLeaves = 0;
	for(vector<PTUNodePtr>::const_iterator nodeIt = id2node.begin(); nodeIt != id2node.end(); ++nodeIt)
		if((*nodeIt)->isLeaf())
			nLeaves++;
	return nLeaves;
}

inline void PTUnrooted::evaluate(const PTUNodePtr& node, const DNASubModel& model) {
	for(int j = 0; j < csLen; ++j) {
//		infoLog << "Evaluating at site " << j << "\r";
		evaluate(node, j, model);
	}
}

inline void PTUnrooted::evaluate(const PTUNodePtr& node, int j, const DNASubModel& model) {
	for(vector<PTUNodePtr>::iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) {
		if(isChild(*child, node) && !isEvaluated(*child, node, j))
			evaluate(*child, node, j, model);
	}
}


inline ostream& PTUnrooted::writeTree(ostream& out, string format) const {
	StringUtils::toLower(format);
	if(format == "newick")
		return writeTreeNewick(out, root) << ";";
	else {
		errorLog << "Cannot write PTUnrooted tree, unknown tree format " << format << std::endl;
		out.setstate(std::ios_base::failbit);
		return out;
	}
}

inline double PTUnrooted::getBranchLength(const PTUNodePtr& u, const PTUNodePtr& v) const {
	return node2length.at(u).at(v);
//	BranchLenMap::const_iterator resultOuter = node2length.find(u);
//	if(resultOuter == node2length.end())
//		return -1;
//	else {
//		boost::unordered_map<PTUNodePtr, double>::const_iterator resultInner = resultOuter->second.find(v);
//		return resultInner == resultOuter->second.end() ? -1 : resultInner->second;
//	}
}

inline Matrix4Xd PTUnrooted::getBranchCost(const PTUNodePtr& u, const PTUNodePtr& v) const {
	CostMap::const_iterator resultOuter = node2cost.find(u);
	if(resultOuter == node2cost.end())
		return Matrix4Xd();
	else {
		boost::unordered_map<PTUNodePtr, Matrix4Xd>::const_iterator resultInner = resultOuter->second.find(v);
		return resultInner == resultOuter->second.end() ? Matrix4Xd() : resultInner->second;
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

inline vector<Matrix4d> PTUnrooted::getModelTransitionSet(string method) const {
	StringUtils::toLower(method);
	if(method == "gojobori")
		return getModelTraningSetGojobori();
	else if(method == "goldman")
		return getModelTraningSetGoldman();
	else
		throw invalid_argument("Unknown DNA substitution model training method '" + method + "'");
}

inline bool PTUnrooted::isTip(const PTUNodePtr& node) {
	if(node->isLeaf())
		return false;
	for(vector<PTUNodePtr>::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child)
		if(isChild(*child, node) && !(*child)->isLeaf())
			return false;
	return true;
}

inline PTUnrooted::PTUNodePtr PhyloTreeUnrooted::firstLeaf(PTUNodePtr node) {
	while(!node->isLeaf())
		node = node->firstChild();
	return node;
}

inline PTUnrooted::PTUNodePtr PhyloTreeUnrooted::lastLeaf(PTUNodePtr node) {
	while(!node->isLeaf())
		node = node->lastChild();
	return node;
}

inline PTUnrooted::PTUNodePtr PhyloTreeUnrooted::randomLeaf(PTUNodePtr node) {
	while(!node->isLeaf()) {
		const vector<PTUNodePtr>& children = node->getChildren();
		node = children[rand() % children.size()];
	}
	return node;
}



} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREEUNROOTED_H_ */
