/*
 * PhyloTreeUnrooted.h
 *  An Unrooted Phylogenic Tree (PTUnrooted)
 *  A PTUnrooted can be evaluated from any node as its root and yields same loglik
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
#include <sstream>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/unordered_map.hpp>

#include "AlphabetFactory.h"
#include "HmmUFOtuConst.h"
#include "ProgLog.h"
#include "StringUtils.h"
#include "DigitalSeq.h"
#include "NewickTree.h"
#include "MSA.h"
#include "DNASubModel.h"
#include "DiscreteGammaModel.h"

namespace EGriceLab {
using std::string;
using std::vector;
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

	class PhyloTreeUnrootedNode;
	typedef PTUnrooted::PhyloTreeUnrootedNode PTUNode;

	class PhyloTreeUnrootedBranch;
	typedef PTUnrooted::PhyloTreeUnrootedBranch PTUBranch;

	typedef shared_ptr<PTUNode> PTUNodePtr; /* use boost shared_ptr to hold node pointers */
	typedef shared_ptr<const PTUNode> PTUNodeConstPtr; /* use boost shared_ptr to hold node pointers */

	typedef shared_ptr<DNASubModel> ModelPtr; /* use boost shared_ptr to hold DNA Sub Model */
	typedef shared_ptr<DiscreteGammaModel> DGammaPtr; /* use boost shared_ptr to hold DiscreteGammapModel */

	typedef boost::unordered_map<PTUNodePtr, boost::unordered_map<PTUNodePtr, PTUBranch> > BranchMap;

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
		 * Construct a PTUNode with a given id, name, sequence, and optionally annotation and annotation-dist
		 */
		PhyloTreeUnrootedNode(long id, const string& name, const DigitalSeq& seq,
				const string& anno = "", double annoDist = 0)
		: id(id), name(name), seq(seq), anno(anno), annoDist(annoDist)
		{ }

		/**
		 * Construct a PTUNode with a given id, name, unobserved seq with given length, and optionally annotation and annotation-dist
		 */
		PhyloTreeUnrootedNode(long id, const string& name, size_t length,
				const string& anno = "", double annoDist = 0)
		: id(id), name(name), seq(AlphabetFactory::getAlphabetByName("DNA"), name),
		  anno(anno), annoDist(annoDist)
		{
			seq.append(length, DegenAlphabet::GAP_SYM);
		}

		/* Member methods */
		/* Getters and Setters */
		const string& getAnno() const {
			return anno;
		}

		double getAnnoDist() const {
			return annoDist;
		}

		/**
		 * Get node annotation, with an optional "Other_dist" posted if too far from its annotation source
		 */
		string getAnnotation(double maxDist = 0) const;

		long getId() const {
			return id;
		}

		const string& getName() const {
			return name;
		}

		string getLabel() const;

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

		/**
		 * load data from a binary input to this node
		 */
		istream& load(istream& in);

		/**
		 * save this node to a binrary output, ignore its edges
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

	class PhyloTreeUnrootedBranch {
		friend class PhyloTreeUnrooted;

	public:
		/** default constructor */
		PhyloTreeUnrootedBranch() : length(0) { }

		/** construct a branch with given length */
		PhyloTreeUnrootedBranch(double length) : length(length) { }

		/** construct a branch with given length and loglik */
		PhyloTreeUnrootedBranch(double length, const Matrix4Xd& loglik) :
			length(length), loglik(loglik)
		{ }

		/** save this branch to a binary output */
		ostream& save(ostream& out) const;

		/** load data from a binary input to this branch */
		istream& load(istream& in);

	private:
		double length; /* branch length */
		Matrix4Xd loglik; /* outgoing message (loglik) of this branch, before convoluting into branch length */
	};

	/* constructors */
	/** Default constructor, do nothing */
	PhyloTreeUnrooted() : csLen(0) {  }

	/** Construct a PTUnrooted from a Newick Tree */
	PhyloTreeUnrooted(const NewickTree& ntree);

public:

	/* member methods */

	/** Get the number of nodes of this tree */
	size_t numNodes() const {
		return id2node.size();
	}

	size_t numEdges() const;

	/** Get number of leaves in this tree */
	size_t numLeaves() const;

	/** Get number of branches in this tree */
	size_t numBranches() const {
		return numLeaves() / 2;
	}

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

	/** add a new edge u<->v to this tree */
	void addEdge(const PTUNodePtr& u, const PTUNodePtr& v) {
		u->neighbors.push_back(v);
		v->neighbors.push_back(u);
	}

	/** remove an edge u<->v to this tree */
	void removeEdge(const PTUNodePtr& u, const PTUNodePtr& v) {
		u->neighbors.erase(std::find(u->neighbors.begin(), u->neighbors.end(), v));
		v->neighbors.erase(std::find(v->neighbors.begin(), v->neighbors.end(), u));
	}

	/**
	 * get branch from u-> v
	 * @throw  out_of_range exception if not exists
	 */
	const PTUBranch& getBranch(const PTUNodePtr& u, const PTUNodePtr& v) const {
		return node2branch.at(u).at(v);
	}

	/**
	 * set branch from u-> v
	 */
	void setBranch(const PTUNodePtr& u, const PTUNodePtr& v, const PTUBranch& w) {
		node2branch[u][v] = w;
	}

	/**
	 * remove the branch from u->v
	 * @return  the old branch
	 */
	void removeBranch(const PTUNodePtr& u, const PTUNodePtr& v) {
		node2branch[u].erase(node2branch[u].find(v));
	}

	/**
	 * get branch length from u -> v
	 * @return  branch length u->v
	 * @throw  out_of_range exception if branch not exists
	 */
	double getBranchLength(const PTUNodePtr& u, const PTUNodePtr& v) const {
		return node2branch.at(u).at(v).length;
	}

	/**
	 * set branch length from u <-> v
	 */
	void setBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, double w) {
		node2branch[u][v].length = node2branch[v][u].length = w;
	}

	/**
	 * get branch loglik of u->v at site j
	 */
	Vector4d getBranchLoglik(const PTUNodePtr& u, const PTUNodePtr& v, int j) const {
		return node2branch.at(u).at(v).loglik.col(j);
	}

	/**
	 * get branch loglik of u->v at all sites
	 */
	const Matrix4Xd& getBranchLoglik(const PTUNodePtr& u, const PTUNodePtr& v) const {
		return node2branch.at(u).at(v).loglik;
	}

	/**
	 * set branch loglik of u->v at site j
	 */
	void setBranchLoglik(const PTUNodePtr& u, const PTUNodePtr& v, int j, const Vector4d& loglik) {
		node2branch[u][v].loglik.col(j) = loglik;
	}

	/**
	 * set branch loglik of u->v at all sites
	 */
	void setBranchLoglik(const PTUNodePtr& u, const PTUNodePtr& v, const Matrix4Xd& loglik) {
		node2branch[u][v].loglik = loglik;
	}

	/** Load sequences from MSA into this this */
	size_t loadMSA(const MSA& msa);

	/** Load tab-delimited annotation file of tree nodes into this tree */
	istream& loadAnnotation(istream& in);

	/** format node names to exclude white spaces and unprintable characters */
	void formatName();

	/** format node annotations to exclude white spaces and unprintable characters */
	void formatAnnotation();

	/**
	 * annotate every node of this tree by checking their neighbors' names
	 */
	void annotate();

	/**
	 * annotate a node, either by itself or by a named nearest neighbor
	 */
	void annotate(const PTUNodePtr& node);

	/**
	 * annotate a node by its parent using DFS
	 */
	void annotate(const PTUNodePtr& u, const PTUNodePtr& v, boost::unordered_map<PTUNodePtr, double>& visited);

	/**
	 * Set the underlying DNA Sub Model as a copy of given model
	 */
	void setModel(const DNASubModel& model) {
		this->model.reset(model.clone());
	}

	/**
	 * Set the underlying DNA Sub Model as a copy of this object
	 */
	void setModel(const DNASubModel* model) {
		this->model.reset(model->clone());
	}

	/**
	 * Get the underlying DNA Sub Model
	 */
	const ModelPtr& getModel() const {
		return model;
	}

	/**
	 * Set the underlying DG Model as a copy of given model
	 */
	void setDGModel(const DiscreteGammaModel& dG) {
		this->dG.reset(dG.clone());
	}

	/**
	 * Set the underlying DNA Sub Model as a copy of this object
	 */
	void setDGModel(const DiscreteGammaModel* dG) {
		this->dG.reset(dG->clone());
	}

	/**
	 * Get the underlying Discrete Gamma Model
	 */
	const DGammaPtr& getDGModel() const {
		return dG;
	}

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
	 * test whether the loglik (message) of node u -> v of all site j has been evaluated
	 */
	bool isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v) const;

	/**
	 * test whether the loglik (message) of node u -> v of site j has been evaluated
	 */
	bool isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v, int j) const;

	/**
	 * initiate the cached incoming loglik of between every node u and every neighbor v
	 */
	void initInLoglik();

	/**
	 * initiate the leaf loglik
	 */
	void initLeafLoglik();

	/**
	 * reset the cached loglik of edge u->v
	 */
	void resetLoglik(const PTUNodePtr& u, const PTUNodePtr& v) {
		node2branch[u][v].loglik.setConstant(INVALID_LOGLIK);
	}

	/**
	 * reset the cached loglik of every node
	 */
	void resetLoglik();

	/**
	 * reset the cached leaf loglik
	 */
	void resetLeafLoglik() {
		leafLoglik.setConstant(INVALID_LOGLIK);
	}

	/**
	 * evaluate the conditional loglik of the jth site of a subtree,
	 * rooted at given node, with a given rate factor r
	 * this is the base for all evaluate/loglik methods
	 * @param node  subtree root
	 * @param j  the jth aligned site
	 * @param r  the rate factor at site j
	 * @return  conditional loglik at the jth site
	 */
	Vector4d loglik(const PTUNodePtr& node, int j, double r);

	/**
	 * evaluate the conditional loglik of the jth site of a subtree, rooted at given node
	 * this method will use either fixed rate (r===1) or a DiscreteGammaModel to evaluate the loglik
	 * according to Yang 1994 (b)
	 * @param node  subtree root
	 * @param j  the jth aligned site
	 * @return  conditional loglik at the jth site
	 */
	Vector4d loglik(const PTUNodePtr& node, int j);

	/**
	 * evaluate the log-likelihood (loglik) of the entire tree
	 * @return  loglik matrix of the entire tree
	 */
	Matrix4Xd loglik() {
		return loglik(root);
	}

	/**
	 * evaluate the log-likelihood (loglik) at the jth site of the entire tree
	 * @return  loglik vector at the jth site
	 */
	Vector4d loglik(int j) {
		return loglik(root, j);
	}

	/**
	 * evaluate the conditional loglik of a subtree, rooted at given node
	 * @param node  subtree root
	 * @return  conditional loglik matrix of the subtree
	 */
	Matrix4Xd loglik(const PTUNodePtr& node);

	/**
	 * evaluate the entire tree
	 */
	void evaluate() {
		evaluate(root);
	}

	/**
	 * evaluate the subtree at given node
	 */
	void evaluate(const PTUNodePtr& node) {
		for(int j = 0; j < csLen; ++j)
			evaluate(node, j);
	}

	/**
	 * evaluate every child of this node at the jth site of a subtree, rooted at given node
	 * but does not calculate the loglik of this node itself
	 * @param node  subtree root
	 * @param j  the jth aligned site
	 */
	void evaluate(const PTUNodePtr& node, int j);

	/**
	 * calculate the subtree loglik at j-th site
	 * evaluate the tree if necessary
	 * this is the basis for all other treeLoglik methods
	 */
	double treeLoglik(const PTUNodePtr& node, int j);

	/**
	 * calculate the loglik of the subtree in a given range [start, end]
	 */
	double treeLoglik(const PTUNodePtr& node, int start, int end);

	/**
	 * calculate the loglik of the subtree in a whole length
	 */
	double treeLoglik(const PTUNodePtr& node);

	/**
	 * calculate the entire tree loglik at j-th site
	 * evaluate the tree if necessary
	 */
	double treeLoglik(int j);

	/**
	 * calculate the entire tree loglik in a given range [start, end]
	 */
	double treeLoglik(int start, int end);

	/**
	 * calculate the entire tree loglik in the whole length
	 */
	double treeLoglik();

	/**
	 * infer the ancestor (or real if a leaf) state (base) of given node and site
	 * @param node  node to infer
	 * @param j  alignment site
	 * @return  the actual observed state if a leaf node,
	 * or inferred state my maximazing the conditional liklihood
	 */
	int inferState(const PTUNodePtr& node, int j);

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

	vector<PTUNodePtr> getLeafHits(const DigitalSeq& seq, double maxPDist,
			int start, int end) const;

	vector<PTUNodePtr> getLeafHits(const DigitalSeq& seq, double maxPDist) const {
		return getLeafHits(seq, maxPDist, 0, csLen - 1);
	}

	/**
	 * get estimated base frequency (pi) using this tree
	 */
	Vector4d getModelFreqEst() const;

	/**
	 * estimate the total number of mutations at given site,
	 * using ML estimation based on conditional liklihoods
	 * @param j  site to estimate
	 * @return  total estimated mutations at this site
	 */
	size_t estimateNumMutations(int j);

	/**
	 * make a copy of subtree with only two nodes and a branch u and v
	 * edges u->v and v->u should has already been evaluated
	 * @return  a new PhyloTreeUnrooted with only two nodes and a branch u->v,
	 * with root set as v
	 */
	PTUnrooted copySubTree(const PTUNodePtr& u, const PTUNodePtr& v) const;

	/**
	 * estimate branch length by comparing the two direction loglik
	 * in given region [start-end]
	 * return the estimated branch length
	 */
	double estimateBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, int start, int end);

	/**
	 * estimate branch length by comparing the two direction loglik
	 * in the entire region
	 * return the estimated branch length
	 */
	double estimateBranchLength(const PTUNodePtr& u, const PTUNodePtr& v) {
		return estimateBranchLength(u, v, 0, csLen - 1);
	}

	/**
	 * iteratively optimize the length of branch u->v using Felsenstein's algorithm
	 * in given CSRegion [start-end]
	 * return the updated branch length v
	 */
	double optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, int start, int end);

	/**
	 * iteratively optimize the length of branch u->v using Felsenstein's algorithm
	 * return the updated branch length v
	 */
	double optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v) {
		return optimizeBranchLength(u, v, 0, csLen - 1);
	}

	/**
	 * place an additional seq (n) at given branch with given initial branch length
	 * by introducing a new internal root r, initially placed at the mid-point between u->v
	 * during placement, tree branch lengths are optimized in the order of
	 * n->r, then u->r/r->v jointly
	 * @param  new seq to be placed
	 * @param u  branch start (u->v)
	 * @param v  branch end (u->v)
	 * @param start  seq start position (non-gap start)
	 * @param end  seq end position (non-gap end)
	 * @param d0  estimiated initial branch length n->r
	 * @return  the modified tree, with r and n appended at the end of all nodes
	 */
	PTUnrooted& placeSeq(const DigitalSeq& seq, const PTUNodePtr& u, const PTUNodePtr& v,
			int start, int end);

	/**
	 * place an additional seq (n) at given branch with given initial branch length
	 * after placement, tree will be re-rooted to the new internal root r,
	 * and new branch length n->r will be optimized according the entire region [0, csLen-1]
	 * @param  new seq to be placed
	 * @param u  branch start (u->v)
	 * @param v  branch end (u->v)
	 * @param d0  estimiated initial branch length n->r
	 * @return  the modified tree, with r and n appended at the end of all nodes
	 */
	PTUnrooted& placeSeq(const DigitalSeq& seq, const PTUNodePtr& u, const PTUNodePtr& v) {
		return placeSeq(seq, u, v, 0, csLen -1);
	}

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
	 * load leaf loglik from a binary input
	 */
	istream& loadLeafLoglik(istream& in);

	/**
	 * save leaf loglik to a binary outout
	 */
	ostream& saveLeafLoglik(ostream& out) const;

	/**
	 * load root information from a binary input
	 */
	istream& loadRoot(istream& in);

	/**
	 * save root information to a binary output
	 */
	ostream& saveRoot(ostream& out) const;

	/**
	 * load root loglik for every node
	 */
	istream& loadRootLoglik(istream& in);

	/**
	 * save root loglik for every node
	 */
	ostream& saveRootLoglik(ostream& out) const;

	/**
	 * load DNA model from a text input
	 */
	istream& loadModel(istream& in);

	/**
	 * save DNA model to a text output
	 */
	ostream& saveModel(ostream& out) const;

	/**
	 * load DiscreteGamma model from a binary input, if any
	 */
	istream& loadDGModel(istream& in);

	/**
	 * save DiscreteGamma model to a binary output, if not NULL
	 */
	ostream& saveDGModel(ostream& out) const;

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

	/* return dot product between a Matrix and a vector, scale the vector if necessary */
	static Vector4d dot_product_scaled(const Matrix4d& X, const Vector4d& V);

	/* return dot product between a pi vector and a loglik vector, scale the second vector if necessary */
	static double dot_product_scaled(const Vector4d& P, const Vector4d& V);

	/* return dot product between two loglik vectors, scale both if necessary */
	static double dot_product_double_scaled(const Vector4d& V1, const Vector4d& V2);

	/* return the rowwise mean of a given matrix at exponential scale, scale each row if neccessary */
	static Vector4d row_mean_exp_scaled(const Matrix4Xd& X);

	/**
	 * format taxonomy name, removes white spaces and unnecessary unnamed taxa prefix
	 * @param taxa  taxa name to be formated
	 * @return  the formated name
	 */
	static string& formatTaxaName(string& taxa);

	static bool isCanonicalName(const string& taxa);

	/* member fields */
private:
	int csLen; /* number of aligned sites */

	PTUNodePtr root; /* root node of this tree */
	vector<PTUNodePtr> id2node; /* indexed tree nodes */

	BranchMap node2branch; /* branch length index storing edge length */
	Matrix4Xd leafLoglik; /* cached 4 X 5 leaf loglik matrix,
						with each column the pre-computed loglik of observing A, C, G, T or - at any given site */

	ModelPtr model; /* DNA Model used to evaluate this tree, needed to be stored with this tree */
	DGammaPtr dG; /* DiscreteGammaModel used to conpensate rate-heterogeinity between alignment sites */

public:
	/* static fields */
	static const double MIN_LOGLIK_EXP;
	static const double INVALID_LOGLIK;

	static const double LOGLIK_REL_EPS;
	static const double BRANCH_EPS;
	static const char ANNO_FIELD_SEP = '\t';
	static const string KINDOM_PREFIX;
	static const string PHYLUM_PREFIX;
	static const string CLASS_PREFIX;
	static const string ORDER_PREFIX;
	static const string FAMILY_PREFIX;
	static const string GENUS_PREFIX;
	static const string SPECIES_PREFIX;
};

inline size_t PTUnrooted::numEdges() const {
	size_t N = 0;
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		N += (*node)->numNeighbors();
	return N;
}

inline size_t PTUnrooted::numLeaves() const {
	size_t N = 0;
	for(vector<PTUNodePtr>::const_iterator nodeIt = id2node.begin(); nodeIt != id2node.end(); ++nodeIt)
		if((*nodeIt)->isLeaf())
			N++;
	return N;
}

inline std::string PTUnrooted::PTUNode::getLabel() const {
	string label;
	std::ostringstream os(label);
	os << anno << annoDist << ';';
	return label;
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

inline bool PTUnrooted::isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v) const {
	BranchMap::const_iterator outerResult = node2branch.find(u);
	if(outerResult != node2branch.end()) {
		unordered_map<PTUNodePtr, PTUBranch>::const_iterator innerResult = outerResult->second.find(v);
		if(innerResult != outerResult->second.end())
			return innerResult->second.loglik.cols() == csLen && /* Matrix is initiated */
					(innerResult->second.loglik.array() != INVALID_LOGLIK).all(); /* values are all valid */
	}
	return false;
}

inline bool PTUnrooted::isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v, int j) const {
	BranchMap::const_iterator outerResult = node2branch.find(u);
	if(outerResult != node2branch.end()) {
		unordered_map<PTUNodePtr, PTUBranch>::const_iterator innerResult = outerResult->second.find(v);
		if(innerResult != outerResult->second.end())
			return innerResult->second.loglik.cols() == csLen && /* Matrix is initiated */
					(innerResult->second.loglik.col(j).array() != INVALID_LOGLIK).all(); /* values are not invalid */
	}
	return false;
}

inline double PTUnrooted::treeLoglik(int j) {
	return treeLoglik(root, j);
}

inline double PTUnrooted::treeLoglik(int start, int end) {
	return treeLoglik(root, start, end);
}

inline double PTUnrooted::treeLoglik() {
	return treeLoglik(root);
}

inline double PTUnrooted::treeLoglik(const PTUNodePtr& node, int j) {
	return dot_product_scaled(model->getPi(), loglik(node, j));
}

inline double PTUnrooted::treeLoglik(const PTUNodePtr& node) {
	return treeLoglik(node, 0, csLen - 1);
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

inline Vector4d PTUnrooted::dot_product_scaled(const Matrix4d& X, const Vector4d& V) {
	Vector4d Y;
	double maxV = V.maxCoeff();
	double scale = maxV != infV && maxV < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV : 0;

	for(Vector4d::Index i = 0; i < Y.rows(); ++i)
		Y(i) = ::log(X.row(i).dot((V.array() + scale).exp().matrix())) - scale;
	return Y;
}

inline double PTUnrooted::dot_product_scaled(const Vector4d& P, const Vector4d& V) {
	double maxV = V.maxCoeff();
	double scale = maxV != inf && maxV < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV : 0;

	return ::log(P.dot((V.array() + scale).exp().matrix())) - scale;
}

inline double PTUnrooted::dot_product_double_scaled(const Vector4d& V1, const Vector4d& V2) {
	double maxV1 = V1.maxCoeff();
	double maxV2 = V2.maxCoeff();
	double scale1 = maxV1 != inf && maxV1 < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV1 : 0;
	double scale2 = maxV2 != inf && maxV2 < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV2 : 0;

	return ::log((V1.array() + scale1).exp().matrix().dot((V2.array() + scale2).exp().matrix())) - scale1 - scale2;
}

inline Vector4d PTUnrooted::row_mean_exp_scaled(const Matrix4Xd& X) {
	/* determine rowwise scaling factors */
	Vector4d scale;
	for(Matrix4Xd::Index i = 0; i < X.rows(); ++i) {
		double maxV = X.row(i).maxCoeff();
		scale(i) = maxV != inf && maxV < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV : 0;
	}
	return (X.colwise() + scale).array().exp().rowwise().mean().log().matrix() - scale;
}

inline void PTUnrooted::formatName() {
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		formatTaxaName((*node)->name);
}

inline void PTUnrooted::formatAnnotation() {
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		formatTaxaName((*node)->anno);
}

inline bool PhyloTreeUnrooted::isCanonicalName(const string& taxa) {
	return !taxa.empty() &&
			(StringUtils::startsWith(taxa, KINDOM_PREFIX) ||
			StringUtils::startsWith(taxa, PHYLUM_PREFIX) ||
			StringUtils::startsWith(taxa, CLASS_PREFIX) ||
			StringUtils::startsWith(taxa, ORDER_PREFIX) ||
			StringUtils::startsWith(taxa, FAMILY_PREFIX) ||
			StringUtils::startsWith(taxa, GENUS_PREFIX) ||
			StringUtils::startsWith(taxa, SPECIES_PREFIX));
}

inline string PTUnrooted::PTUNode::getAnnotation(double maxDist) const {
	if(annoDist <= maxDist)
		return anno;
	char dist[32]; /* _ and numbers */
	sprintf(dist, ";Other_%f", annoDist);
	return anno + dist;
}

} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREEUNROOTED_H_ */
