/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
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
#include <set>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <cassert>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/iterator.hpp>
#include <boost/lexical_cast.hpp>

#include "AlphabetFactory.h"
#include "HmmUFOtuConst.h"
#include "ProgLog.h"
#include "StringUtils.h"
#include "DigitalSeq.h"
#include "NewickTree.h"
#include "MSA.h"
#include "DNASubModel.h"
#include "DiscreteGammaModel.h"
#include "EGMath.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::string;
using std::vector;
using std::set;
using std::istream;
using std::ostream;
using Eigen::Matrix4Xd;
using Eigen::Matrix4d;
using Eigen::RowVectorXd;
using boost::shared_ptr;
using boost::unordered_map;
using boost::unordered_set;

class PhyloTreeUnrooted; /* forward declaration */

typedef PhyloTreeUnrooted PTUnrooted;

class PhyloTreeUnrooted {
public:
	/* nested types and enums */
	/** taxon levels */
	enum TaxonLevel {
		/* Domain, */ Kindom, Phylum, Class, Order, Family, Genus, Species
	};

	/** prior probability types */
	enum PRIOR_TYPE {
		UNIFORM,
		HEIGHT
	};

	class PhyloTreeUnrootedNode;
	typedef PTUnrooted::PhyloTreeUnrootedNode PTUNode;

	class PhyloTreeUnrootedBranch;
	typedef PTUnrooted::PhyloTreeUnrootedBranch PTUBranch;

	class PTPlacement;

	typedef shared_ptr<PTUNode> PTUNodePtr; /* use boost shared_ptr to hold node pointers */
	typedef shared_ptr<const PTUNode> PTUNodeConstPtr; /* use boost shared_ptr to hold node pointers */

	typedef shared_ptr<DNASubModel> ModelPtr; /* use boost shared_ptr to hold DNA Sub Model */
	typedef shared_ptr<DiscreteGammaModel> DGammaPtr; /* use boost shared_ptr to hold DiscreteGammapModel */

	typedef boost::unordered_map<PTUNodePtr, boost::unordered_map<PTUNodePtr, PTUBranch> > BranchMap;
	typedef boost::unordered_map<PTUNodePtr, double> HeightMap;

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
		explicit PhyloTreeUnrootedNode(long id, const string& name)
		: id(id), name(name), annoDist(0) {  }

		/**
		 * Construct a PTUNode with a given id, name, annotation and annotation-dist
		 */
		PhyloTreeUnrootedNode(long id, const string& name,
				const string& anno, double annoDist)
		: id(id), name(name), anno(anno), annoDist(annoDist)
		{ }

		/**
		 * Construct a PTUNode with a given id, name and sequence
		 */
		PhyloTreeUnrootedNode(long id, const string& name, const DigitalSeq& seq)
		: id(id), name(name), seq(seq), annoDist(0)
		{ }

		/**
		 * Construct a PTUNode with a given id, name, sequence, annotation and annotation-dist
		 */
		PhyloTreeUnrootedNode(long id, const string& name, const DigitalSeq& seq,
				const string& anno, double annoDist)
		: id(id), name(name), seq(seq), anno(anno), annoDist(annoDist)
		{ }

		/* Member methods */
		/* Getters and Setters */
		const string& getAnno() const {
			return anno;
		}

		/**
		 * Get node taxon annotation, with an optional "Other" suffix if too far from its annotation source
		 */
		string getTaxon(double maxDist = inf) const;

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

		double getAnnoDist() const {
			return annoDist;
		}

		void setAnno(const string& anno) {
			this->anno = anno;
		}

		void setAnnoDist(double annoDist) {
			this->annoDist = annoDist;
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
			return parent == nullNode;
		}

		/** test whether this node is parent of another node */
		bool isParent(const PTUNodePtr& other) const {
			return other != nullNode && this == other->parent.get();
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
		 * return nullNode if not exists
		 */
		PTUNodePtr firstChild() const {
			for(vector<PTUNodePtr>::const_iterator child = neighbors.begin(); child != neighbors.end(); ++child)
				if(isParent(*child)) // this is really a child
					return *child;
			return nullNode;
		}

		/**
		 * get last child of this node
		 * return nullNode if not exists
		 */
		PTUNodePtr lastChild() const {
			for(vector<PTUNodePtr>::const_reverse_iterator child = neighbors.rbegin(); child != neighbors.rend(); ++child)
				if(isParent(*child)) // this is really a child
					return *child;
			return nullNode;
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
		 * save this node to a binary output, ignore its edges
		 */
		ostream& save(ostream& out) const;

	private:
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
		PhyloTreeUnrootedBranch() {  }

		/** construct a branch with given length */
		explicit PhyloTreeUnrootedBranch(double length) : length(length) {  }

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

	/**
	 * A simple POD type for store PT seed locations
	 */
	struct PTLoc {
		/* constructors */
		/** construct from given data */
		PTLoc(int start, int end, long id, double dist)
		: start(start), end(end), id(id), dist(dist)
		{  }

		/* non-member friend functions */
		friend bool operator<(const PTLoc& lhs, const PTLoc& rhs);

		/* member fields */
		int start; /* 0-based aln start */
		int end;   /* 0-based aln end */
		long id;   /* node ID */
		double dist; /* p-Dist to node */
	};

	/**
	 * A candidate Phylogenetic Tree Placement to store placement information
	 */
	struct PTPlacement {
		/* constructors */
	//	/** default constructor */
		PTPlacement() : start(0), end(0), cNode(NULL), pNode(NULL), aNode(NULL),
				wuv(nan), ratio(nan), wnr(nan), loglik(nan), height(nan), qPlace(nan), qTaxon(nan)
		{  }

		/** construct a placement with basic info and optionally auxilary info */
		PTPlacement(int start, int end,
				const PTUnrooted::PTUNodePtr& cNode, const PTUnrooted::PTUNodePtr& pNode, const PTUnrooted::PTUNodePtr& aNode,
				double wuv, double ratio, double wnr, double loglik,
				double height = 0, double qPlace = 0, double qTaxonomy = 0)
		: start(start), end(end), cNode(cNode), pNode(pNode), aNode(aNode),
		  wuv(wuv), ratio(ratio), wnr(wnr), loglik(loglik), height(height), qPlace(qPlace), qTaxon(qTaxonomy)
		{  }

		/** destructor */
		virtual ~PTPlacement() {  }

		/** member methods */
		long getTaxonId() const {
			if(aNode != NULL)
				return aNode->getId();
			else
				return UNASSIGNED_TAXONID;
		}

		string getTaxonName() const {
			if(aNode != NULL)
				return aNode->getAnno();
			else
				return UNASSIGNED_TAXONNAME;
		}

		string getId() const {
			if(cNode != NULL && pNode != NULL)
				return boost::lexical_cast<string> (cNode->getId()) + "->" + boost::lexical_cast<string> (pNode->getId());
			else
				return UNASSIGNED_ID;
		}

		bool isValidPlace() const {
			return isParent(pNode, cNode);
		}

		/** calculate prior probability of a placement given a prior type in log-scale */
		double logPriorPr(PRIOR_TYPE type) const;

		/** calculate prior proability of a placement given a prior type */
		double priorPr(PRIOR_TYPE type) const {
			return ::exp(logPriorPr(type));
		}

		/** get the annotation distance of this placement */
		double getAnnoDist() const {
			assert(isValidPlace());
			return aNode == cNode ? wuv * ratio + wnr : (1 - ratio) * wuv + wnr;
		}

		/** write this PTPlacement to text output */
		ostream& write(ostream& out) const;

		/** get segment tree loglik at given region */
//		double segLoglik(int start, int end) const {
//			return treeLoglik.segment(start, end - start + 1).sum();
//		}

		/* non-member functions */
		friend bool compareByLoglik(const PTPlacement& lhs, const PTPlacement& rhs);
		friend bool compareByQTaxon(const PTPlacement& lhs, const PTPlacement& rhs);
		friend bool compareByQPlace(const PTPlacement& lhs, const PTPlacement& rhs);
		friend ostream& operator<<(ostream& out, const PTPlacement& place);

		/** member fields */
		int start; /* 0-based align start */
		int end;   /* 0-based align end */
		PTUnrooted::PTUNodePtr cNode;  /* child node */
		PTUnrooted::PTUNodePtr pNode;  /* parent node */
		PTUnrooted::PTUNodePtr aNode;  /* assigned/annotated node, can be either cNode or pNode */
		double wuv;   /* placement branch length */
		double ratio; /* placement ratio */
		double wnr;   /* new branch length */
		double loglik;
//		double annoDist;
		double height;
		double qPlace;
		double qTaxon;
//		VectorXd treeLoglik; /* optional entire placement tree loglik at every site */

		/** static member fields */
		static const int MAX_Q = 250; /* maximum allowed Q value */
		static const long UNASSIGNED_TAXONID = -1;
		static const string UNASSIGNED_TAXONNAME;
		static const double UNASSIGNED_LOGLIK;
		static const string UNASSIGNED_ID;
		static const double UNASSIGNED_POSTQ;
		static const double UNASSIGNED_DIST;
		static const double UNASSIGNED_RATIO;

		static const string TSV_HEADER;
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
		return numNodes();
	}

	/** get number of aligned sites */
	int numAlignSites() const {
		return csLen;
	}

	/** get root node */
	const PTUNodePtr& getRoot() const {
		return root;
	}

	/** get MSA2Node index */
	const map<unsigned, PTUNodePtr>& getMSA2NodeIndex() const {
		return msaId2node;
	}

	/** get Node2MSA index */
	const map<PTUNodePtr, unsigned>& getNode2MSAIndex() const {
		return node2msaId;
	}

	/** get node by MSA id */
	PTUNodePtr getNodeByMSAId(unsigned id) const {
		return msaId2node.at(id);
	}

	/** get MSA id by node */
	unsigned getMSAIdByNode(const PTUNodePtr& node) const {
		return node2msaId.at(node);
	}

	/** get all nodes */
	std::vector<PTUNodePtr> getNodes() const {
		return id2node;
	}

	/** get node i */
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
	 * add edgeID from given pair of branches
	 * @return  the childId if they are parent/child, or -1 if not
	 */
	long getEdgeID(const PTUNodePtr& u, const PTUNodePtr& v) const {
		if(isChild(u, v))
			return u->id;
		else if(isChild(v, u))
			return v->id;
		else
			return -1;
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

	/**
	 * get node height given node ptr
	 */
	double getHeight(const PTUNodePtr& node) const {
		return node2height.at(node);
	}

	/**
	 * get node height given node id
	 */
	double getHeight(long id) const {
		return getHeight(id2node[id]);
	}

	/**
	 * get all node heights
	 */
	 const HeightMap& getHeights() const {
		 return node2height;
	 }

	/** Load sequences from MSA into this tree
	 * @param msa  MSA data to load
	 * @return  number of loaded nodes, or -1 if error happend
	 */
	unsigned loadMSA(const MSA& msa);

	/** Load tab-delimited annotation file of tree nodes into this tree */
	istream& loadAnnotation(istream& in);

	/** format node names to exclude white spaces and unprintable characters */
	void formatName();

	/** format node annotations to exclude white spaces and unprintable characters */
	void formatAnnotation();

	/**
	 * annotate every node of this tree
	 */
	void annotate(const string& rootName = DEFAULT_ROOT_NAME);

	/**
	 * annotate a node, either by itself or by a named nearest neighbor
	 */
	void annotate(const PTUNodePtr& node, const string& rootName = DEFAULT_ROOT_NAME);

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
	 * test whether among-site varation is enabled
	 */
	bool isVar() const {
		return dG != nulldG;
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
	 * save PTUnrooted to binary output
	 */
	ostream& save(ostream& out) const;

	/** load PTUnrooted from a binary input */
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

	/** calculate all node height at current root */
	void calcNodeHeight();

	/**
	 * fix tree branch length with non-positive values
	 * @param minLen  minimum branch length to use, must be non-negative
	 */
	void fixBranchLength(double minLen = BRANCH_EPS);

	/**
	 * test whether the loglik (message) of node u -> v of all site j has been evaluated
	 */
	bool isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v) const;

	/**
	 * test whether the loglik (message) of node u -> v of site j has been evaluated
	 */
	bool isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v, int j) const;

	/**
	 * test whether the loglik of node u->v of this region is evaluated
	 */
	bool isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v, int start, int end) const;

	/**
	 * initiate the cached incoming loglik of between every node u and every neighbor v
	 */
	void initBranchLoglik();

	/**
	 * initiate the cached root loglik
	 */
	void initRootLoglik() {
		node2branch[root][nullNode].loglik = Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK);
	}

	/**
	 * update the cached root loglik
	 */
	void updateRootLoglik();

	/**
	 * reset the cached loglik of edge u->v
	 */
	void resetLoglik(const PTUNodePtr& u, const PTUNodePtr& v) {
		node2branch[u][v].loglik.setConstant(INVALID_LOGLIK);
	}

	/**
	 * reset the cached loglik of edge u->v at given region
	 */
	void resetLoglik(const PTUNodePtr& u, const PTUNodePtr& v, int start, int end) {
		node2branch[u][v].loglik.block(0, start, 4, end - start + 1).setConstant(INVALID_LOGLIK);
	}

	/**
	 * reset the cached loglik of every node
	 */
	void resetBranchLoglik();

	/**
	 * reset the cached root loglik
	 */
	void resetRootLoglik() {
		node2branch[root][nullNode].loglik.setConstant(INVALID_LOGLIK);
	}

	/**
	 * evaluate the convoluted conditional loglik of the jth site of a subtree,
	 * rooted at given node, with a given rate factor r
	 * this is the base for all evaluate/loglik methods
	 * @param node  subtree root
	 * @param j  the jth aligned site
	 * @param r  the rate factor at site j
	 * @return  convoluted conditional loglik at the jth site
	 */
	Vector4d loglikConv(const PTUNodePtr& node, int j, double r = 1) const;

	/**
	 * evaluate the convoluted conditional loglik of the jth site of a subtree,
	 * rooted at given node, with a given rate factor r
	 * this is the base for all evaluate/loglik methods
	 * @param node  subtree root
	 * @param j  the jth aligned site
	 * @param r  the rate factor at site j
	 * @return  convoluted conditional loglik at the jth site
	 */
	Vector4d loglik(const PTUNodePtr& node, int j) const;

	/**
	 * evaluate the log-likelihood (loglik) of the entire tree
	 * @return  loglik matrix of the entire tree
	 */
	Matrix4Xd loglik() const {
		return loglik(root);
	}

	/**
	 * evaluate the log-likelihood (loglik) at the jth site of the entire tree
	 * @return  loglik vector at the jth site
	 */
	Vector4d loglik(int j) const {
		return loglik(root, j);
	}

	/**
	 * evaluate the conditional loglik of a subtree, rooted at given node
	 * @param node  subtree root
	 * @return  conditional loglik matrix of the subtree
	 */
	Matrix4Xd loglik(const PTUNodePtr& node) const;

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
		evaluate(node, 0, csLen - 1);
	}

	/**
	 * evaluate the subtree at given node at given region
	 */
	void evaluate(const PTUNodePtr& node, int start, int end);

	/**
	 * calculate the loglike of the subtree at site j
	 */
	double treeLoglik(const PTUNodePtr& node, int j) const {
		return treeLoglik(model->getPi(), getBranchLoglik(node, node->parent), j);
	}

	/**
	 * calculate the loglik of the subtree in a given range [start, end]
	 */
	double treeLoglik(const PTUNodePtr& node, int start, int end) const;

	/**
	 * calculate the loglik of the subtree in a whole length
	 */
	double treeLoglik(const PTUNodePtr& node) const {
		return treeLoglik(node, 0, csLen - 1);
	}

	/**
	 * calculate the tree loglik at given site for root node
	 */
	double treeLoglik(int j) const {
		return treeLoglik(root, j);
	}

	/**
	 * calculate the entire tree loglik in a given range [start, end]
	 */
	double treeLoglik(int start, int end) const {
		return treeLoglik(root, start, end);
	}

	/**
	 * calculate the entire tree loglik in the whole length
	 */
	double treeLoglik() const {
		return treeLoglik(root);
	}

	/**
	 * infer the ancestor (or real if a leaf) state (base) of given node and site
	 * @param node  node to infer
	 * @param j  alignment site
	 * @return  the actual observed state if a leaf node,
	 * or inferred state my maximazing the conditional likelihood
	 */
	int8_t inferState(const PTUNodePtr& node, int j) const {
		return inferState(node, node->parent, j);
	}

	/**
	 * infer the ancestor (or real if a leaf) state (base) of given branch and site
	 * @param u  node to infer
	 * @param v  direction to infer
	 * @param j  alignment site
	 * @return  the actual observed state if a leaf node,
	 * or inferred state my maximazing the conditional likelihood
	 */
	int8_t inferState(const PTUNodePtr& u, const PTUNodePtr& v, int j) const;

	/**
	 * Infer the ancestor sequence of this node
	 * the underlying seq will be resized and modified during inferring
	 * before inferring, the conditional likelihood of this sequence should have been evaluated
	 * it will not modify the seq if it is already inferred or assigned
	 * return true if this node is actually inferred
	 */
	void inferSeq(const PTUNodePtr& node);

	/** Infer all non-leaf node in a tree */
	void inferSeq();

	/** convert this PTUnrooted subtree into NewickTree */
	NewickTree convertToNewickTree(const PTUNodePtr& subtree, const string& prefix = "") const;

	/** convert this PTUUnrooted subtree into NewickTree, only for a subset of nodes */
	NewickTree convertToNewickTree(const PTUNodePtr& subtree,
			const boost::unordered_set<PTUNodePtr>& subset, const string& prefix = "") const;

	/** convert this PTUnrooted tree into NewickTree */
	NewickTree convertToNewickTree(const string& prefix = "") const {
		return convertToNewickTree(root, prefix);
	}

	/** convert this PTUUnrooted tree into NewickTree, only for a subset of nodes */
	NewickTree convertToNewickTree(const boost::unordered_set<PTUNodePtr>& subset, const string& prefix = "") const {
		return convertToNewickTree(root, subset, prefix);
	}

	/**
	 * convert this subtree to JPlace tree string, which has additional edge number information and violate the Newick grammar
	 */
	string toJPlaceTreeStr(const PTUNodePtr& node) const;

public:
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

	/**
	 * estimate the total number of mutations at given site,
	 * using ML estimation based on conditional likelihoods
	 * @param j  site to estimate
	 * @return  total estimated mutations at this site
	 */
	size_t estimateNumMutations(int j) const;

	/** get leaf loglik at site j assuming its seq is the given seq */
	Vector4d getLeafLoglik(const DigitalSeq& seq, int j) const;

	/**
	 * get leaf loglik matrix but only evaluate the value in given region [start, end]
	 * while values outside the region is set to -inf
	 */
	Matrix4Xd getLeafLoglik(const DigitalSeq& seq, int start, int end) const;

	Matrix4Xd getLeafLoglik(const DigitalSeq& seq) const {
		return getLeafLoglik(seq, 0, csLen - 1);
	}

	/**
	 * make a copy of subtree with only two nodes and a branch u and v,
	 * but ignore any assigned sequence
	 * edges u->v and v->u should has already been evaluated
	 * @return  a new PhyloTreeUnrooted with only two nodes and a branch u->v, and their branch loglik
	 * with root set as v
	 */
	PTUnrooted copySubTree(const PTUNodePtr& u, const PTUNodePtr& v) const;

	double estimateBranchLength(const PTUNodePtr& u, const PTUNodePtr& v,
			int start, int end, const string& method = "weighted") const
	{
		return estimateBranchLength(getBranchLoglik(u, v), getBranchLoglik(v, u), start, end, method);
	}

	/**
	 * estimate branch length by comparing the two direction loglik
	 * in the entire region
	 * this method is not responsible for re-evluate the tree after branch-length is modified
	 *
	 * return the estimated branch length
	 */
	double estimateBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, const string& method = "weighted") const {
		return estimateBranchLength(u, v, 0, csLen - 1, method);
	}

	/**
	 * iteratively optimize the length of branch u->v using Felsenstein's algorithm
	 * in given CSRegion [start-end], while the max length is optionally constrained
	 * this method will use the original branch length as its initial guess
	 * return the updated branch length v
	 */
	double optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, int start, int end, double maxL = inf);

	/**
	 * iteratively optimize the length of branch u->v using Felsenstein's algorithm
	 * return the updated branch length v
	 */
	double optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, double maxL = inf) {
		return optimizeBranchLength(u, v, 0, csLen - 1, inf);
	}

	/**
	 * iteratively optimize the branch n->r, u->r and v->r jointly
	 * in given CSRegion [start-end], so the total length wur + wrv won't changed, and wnr update accordingly
	 * before calling this method, all incoming loglik n->r, u->r and v->r should be evaluated
	 * return the optimized branch ratio (wur / wrv)
	 */
	double optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, const PTUNodePtr& r, const PTUNodePtr& n,
			int start, int end);

	/**
	 * iteratively optimize the branch n->r, u->r and v->r jointly
	 * in the entire seq
	 * return the optimized branch ratio (wur / wrv)
	 */
	double optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, const PTUNodePtr& r, const PTUNodePtr& n) {
		return optimizeBranchLength(u, v, r, n, 0, csLen - 1);
	}

	/**
	 * estimate placement given a potential placement loc
	 * the tree breaches will be only evaluated in one path in the order of wnr -> wur -> wvr
	 * and the ratio, wnr and loglik will be estimated
	 * @param  new seq to be estimate placement
	 * @param place  placement holder
	 * @return  modified placement
	 */
	PTPlacement estimateSeq(const DigitalSeq& seq, const PTLoc& loc, const string& method = "weighted") const;

	/**
	 * place an additional seq (n) at given branch in given region [start,end]
	 * by introducing a new internal root r, which will be placed at the initial ratio0 = wur / (wuv)
	 * and the new branch n->r set to initial length wnr0
	 * then all three new branches will be optimized jointly
	 * the modified tree will have r as its new root
	 * @param seq  new seq to be placed
	 * @param u  branch start (u->v)
	 * @param v  branch end (u->v)
	 * @param start  seq start position (non-gap start)
	 * @param end  seq end position (non-gap end)
	 * @param ratio0  insert point
	 * @param wnr0  new branch initial length
	 * @return  the final treeLoglik after placing this read
	 */
	double placeSeq(const DigitalSeq& seq, const PTUNodePtr& u, const PTUNodePtr& v,
			int start, int end, double ratio0, double wnr0);

	/**
	 * place an additional seq (n) at given placement position,
	 * by copying a subtree at given position then do placement,
	 * which will not affect the oroginal tree
	 * after placement, all branch lengths, ratio and loglik will be updated
	 * @param seq  new seq to be placed at a copy of subtree
	 * @param place  given placement position
	 * @param maxHeight  maximum height of the annotation source node
	 * @return  the subtree used for this placement
	 */
	PTUnrooted placeSeq(const DigitalSeq& seq, PTPlacement& place, double maxHeight = inf) const;

	/**
	 * place an additional seq (n) at given branch in the entire seq region
	 * by introducing a new internal root r, which will be placed at the mid-point between u->v
	 * the new branch n->r will be optimized, and direction loglik will be evaluated
	 * @param  new seq to be placed
	 * @param u  branch start (u->v)
	 * @param v  branch end (u->v)
	 * @return  the final treeLoglik after placing this read
	 */
	double placeSeq(const DigitalSeq& seq, const PTUNodePtr& u, const PTUNodePtr& v,
			double ratio0, double wnr0) {
		return placeSeq(seq, u, v, 0, csLen -1, ratio0, wnr0);
	}

	/**
	 * get posterial consensus sequence (CS) of a node using observed count data,
	 * based on Dirichlet Density model and a given prior
	 * @param node  node to infer CS
	 * @param count  observed base frequency matrix for this node
	 * @param alpha  consenstraction prameter of the Dirichlet Distribution as alpha = Sigma(alpha1..K)
	 */
	DigitalSeq inferPostCS(const PTUNodePtr& node, const Matrix4Xd& count, double alpha) const;

	/**
	 * get posterial consensus sequence (CS) of a node using observed count and gap count
	 * based on Dirichlet Density model and a given prior
	 * @param node  node to infer CS
	 * @param count  observed base frequency matrix for this node
	 * @param gap  observed gap frequency for this node
	 * @param alpha  consenstraction prameter of the Dirichlet Distribution as alpha = Sigma(alpha1..K)
	 */
	DigitalSeq inferPostCS(const PTUNodePtr& node, const Matrix4Xd& count, const RowVectorXd& gap, double alpha) const;

private:
	/** save msaId2node index to a binary output */
	ostream& saveMSAIndex(ostream& out) const;

	/** load msaId2node index from a binary input */
	istream& loadMSAIndex(istream& in);

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
	 * load node height from a binary input
	 */
	istream& loadNodeHeight(istream& in);

	/**
	 * save node height to a binary input
	 */
	ostream& saveNodeHeight(ostream& out) const;

	/**
	 * load root information from a binary input
	 */
	istream& loadRoot(istream& in);

	/**
	 * save root information to a binary output
	 */
	ostream& saveRoot(ostream& out) const;

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
		return c != nullNode && c->parent == p;
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

	/*
	 * return dot product between two matrix in given region [start, end],
	 * and leave all other region values unspecified,
	 * scale the second matrix if necessary
	 */
	static Matrix4Xd dot_product_scaled(const Matrix4d& X, const Matrix4Xd& V, int start, int end);

	/* return dot product between two matrix, scale the second matrix if necessary */
	static Matrix4Xd dot_product_scaled(const Matrix4d& X, const Matrix4Xd& V) {
		return dot_product_scaled(X, V, 0, V.cols() - 1);
	}

	/* return dot product between a Matrix and a vector, scale the vector if necessary */
	static Vector4d dot_product_scaled(const Matrix4d& X, const Vector4d& V);

	/* return dot product between a pi vector and a loglik vector, scale the second vector if necessary */
	static double dot_product_scaled(const Vector4d& P, const Vector4d& V);

	/* return dot product between two loglik vectors, scale both if necessary */
	static double dot_product_double_scaled(const Vector4d& V1, const Vector4d& V2);

	/* return the rowwise mean of a given matrix at exponential scale, scale each row if neccessary */
	static Vector4d row_mean_exp_scaled(const Matrix4Xd& X);

	/** get taxon prefix by their level */
	static string taxonLevel2prefix(TaxonLevel level);

	/** test whether this taxon subpart is a carnonical name at any level */
	static bool isCanonicalName(const string& taxon);

	/** test whether this taxon subpart is a carnonical name at a given level */
	static bool isCanonicalName(const string& taxon, TaxonLevel level);

	/** test whether this taxon name is in full canonical format */
	static bool isFullCanonicalName(const string& taxon);

	/** test whether this taxon name is in partial or full canonical format */
	static bool isPartialCanonicalName(const string& taxon);

	/**
	 * format taxonomy name, removes white spaces and unnecessary unnamed taxon prefix
	 * @param taxon  taxon name to be formated
	 * @return  the formated name. which is empty or carnonical like 'k__xxx;p__xxx;c__xxx'
	 */
	static string formatTaxonName(const string& taxon);

	/**
	 * Infer the base based on a given loglik vector
	 */
	static int8_t inferState(const Vector4d& loglik);

	/** Infer the relative weight of each state */
	static Vector4d inferWeight(const Vector4d& loglik);

	/** Estimate branch length using two incoming loglik Matrix in given region [start, end] */
	static double estimateBranchLength(const Matrix4Xd& U, const Matrix4Xd& V,
			int start, int end, const string& method = "weighted");

	/** Estimate branch length using two incoming loglik Matrix, using unweighted difference by ML infeerring */
	static double estimateBranchLengthUnweighted(const Matrix4Xd& U, const Matrix4Xd& V,
			int start, int end);

	/** Estimate branch length using two incoming loglik Matrix, using unweighted difference by ML infeerring */
	static double estimateBranchLengthWeighted(const Matrix4Xd& U, const Matrix4Xd& V,
			int start, int end);

	static double treeLoglik(const Vector4d& pi, const Matrix4Xd& X, int j) {
		return dot_product_scaled(pi, X.col(j));
	}

	static double treeLoglik(const Vector4d& pi, const Matrix4Xd& X, int start, int end);

	static double treeLoglik(const Vector4d& pi, const Matrix4Xd& X) {
		return treeLoglik(pi, X, 0, X.cols() - 1);
	}

	/** initiate the leaf loglik matrix */
	static Matrix4d initLeafMat();

	static boost::unordered_set<PTUNodePtr> getAncestors(const boost::unordered_set<PTUNodePtr>& subset);

	/* member fields */
private:
	int csLen; /* number of aligned sites */

	PTUNodePtr root; /* root node of this tree */
	vector<PTUNodePtr> id2node; /* indexed tree nodes */
	map<unsigned, PTUNodePtr> msaId2node; /* original id in MSA to node map */
	map<PTUNodePtr, unsigned> node2msaId; /* node to original id in MSA map */

	BranchMap node2branch; /* branch length index storing edge length */
	HeightMap node2height; /* node hight (distance to closest leaf */

	ModelPtr model; /* DNA Model used to evaluate this tree, needed to be stored with this tree */
	DGammaPtr dG; /* DiscreteGammaModel used to conpensate rate-heterogeinity between alignment sites */

	static const DGammaPtr nulldG; /* internal null dG model */
	static const PTUNodePtr nullNode; /* internal null node */

public:
	/* static fields */
	static const double MIN_LOGLIK_EXP;
	static const double INVALID_LOGLIK;

	static const double LOGLIK_REL_EPS;
	static const double BRANCH_EPS;
	static const int MAX_ITER = 100;
	static const char ANNO_FIELD_SEP = '\t';
	static const string DOMAIN_PREFIX;
	static const string KINDOM_PREFIX;
	static const string PHYLUM_PREFIX;
	static const string CLASS_PREFIX;
	static const string ORDER_PREFIX;
	static const string FAMILY_PREFIX;
	static const string GENUS_PREFIX;
	static const string SPECIES_PREFIX;

	static const string DEFAULT_ROOT_NAME;

	static const Matrix4d leafMat; /* cached 4 X  leaf loglik matrix,
						with each column the pre-computed loglik of observing A, C, G, T at any given site */
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

inline bool PTUnrooted::isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v) const {
	return node2branch.count(u) != 0 &&
			node2branch.at(u).count(v) != 0 &&
			(node2branch.at(u).at(v).loglik.array() != INVALID_LOGLIK).all();
}

inline bool PTUnrooted::isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v, int j) const {
	return node2branch.count(u) != 0 &&
			node2branch.at(u).count(v) != 0 &&
			(node2branch.at(u).at(v).loglik.col(j).array() != INVALID_LOGLIK).all();
}

inline bool PTUnrooted::isEvaluated(const PTUNodePtr& u, const PTUNodePtr& v, int start, int end) const {
	return node2branch.count(u) != 0 &&
			node2branch.at(u).count(v) != 0 &&
			(node2branch.at(u).at(v).loglik.block(0, start, 4, end - start + 1).array() != INVALID_LOGLIK).all();
}

inline Vector4d PTUnrooted::getLeafLoglik(const DigitalSeq& seq, int j) const {
	int8_t base = seq[j];
	if(base >= 0)
		return leafMat.col(base);
	else
		return model->getPi().array().log();
}

inline Matrix4Xd PTUnrooted::getLeafLoglik(const DigitalSeq& seq, int start, int end) const {
	assert(seq.length() == csLen);
	Matrix4Xd loglik = Matrix4Xd::Constant(4, csLen, infV);
	for(int j = start; j <= end; ++j)
		loglik.col(j) = getLeafLoglik(seq, j);
	return loglik;
}

inline int8_t PhyloTreeUnrooted::inferState(const PTUNodePtr& u, const PTUNodePtr& v, int j) const {
	assert(isParent(v, u) || isParent(u, v));
	return PTUnrooted::inferState(getBranchLoglik(u, v, j));
}

inline void PhyloTreeUnrooted::inferSeq() {
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		if(!(*node)->isLeaf()) /* not a leaf node */
			inferSeq(*node);
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

inline Matrix4Xd PTUnrooted::dot_product_scaled(const Matrix4d& X, const Matrix4Xd& Y, int start, int end) {
	Matrix4Xd Z(4, Y.cols());
	for(Matrix4Xd::Index j = start; j <= end; ++j)
		Z.col(j) = dot_product_scaled(X, static_cast<const Vector4d&> (Y.col(j)));
	return Z;
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
	double scale = maxV != infV && maxV < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV : 0;

	return ::log(P.dot((V.array() + scale).exp().matrix())) - scale;
}

inline double PTUnrooted::dot_product_double_scaled(const Vector4d& V1, const Vector4d& V2) {
	double maxV1 = V1.maxCoeff();
	double maxV2 = V2.maxCoeff();
	double scale1 = maxV1 != infV && maxV1 < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV1 : 0;
	double scale2 = maxV2 != infV && maxV2 < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV2 : 0;

	return ::log((V1.array() + scale1).exp().matrix().dot((V2.array() + scale2).exp().matrix())) - scale1 - scale2;
}

inline Vector4d PTUnrooted::row_mean_exp_scaled(const Matrix4Xd& X) {
	/* determine rowwise scaling factors */
	Vector4d scale;
	for(Matrix4Xd::Index i = 0; i < X.rows(); ++i) {
		double maxV = X.row(i).maxCoeff();
		scale(i) = maxV != infV && maxV < MIN_LOGLIK_EXP ? MIN_LOGLIK_EXP - maxV : 0;
	}
	return (X.colwise() + scale).array().exp().rowwise().mean().log().matrix() - scale;
}

inline string PTUnrooted::taxonLevel2prefix(TaxonLevel level) {
	switch(level) {
//	case Domain:
//		return DOMAIN_PREFIX;
	case Kindom:
		return KINDOM_PREFIX;
	case Phylum:
		return PHYLUM_PREFIX;
	case Class:
		return CLASS_PREFIX;
	case Order:
		return ORDER_PREFIX;
	case Family:
		return FAMILY_PREFIX;
	case Genus:
		return GENUS_PREFIX;
	case Species:
		return SPECIES_PREFIX;
	default:
		return "";
	}
}

inline void PTUnrooted::formatName() {
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		(*node)->name = formatTaxonName((*node)->name);
}

inline void PTUnrooted::formatAnnotation() {
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		(*node)->anno = formatTaxonName((*node)->anno);
}

inline bool PTUnrooted::isCanonicalName(const string& taxon) {
	return  taxon.length() > 3 &&
			(StringUtils::startsWith(taxon, DOMAIN_PREFIX) ||
			StringUtils::startsWith(taxon, KINDOM_PREFIX) ||
			StringUtils::startsWith(taxon, PHYLUM_PREFIX) ||
			StringUtils::startsWith(taxon, CLASS_PREFIX) ||
			StringUtils::startsWith(taxon, ORDER_PREFIX) ||
			StringUtils::startsWith(taxon, FAMILY_PREFIX) ||
			StringUtils::startsWith(taxon, GENUS_PREFIX) ||
			StringUtils::startsWith(taxon, SPECIES_PREFIX));
}

inline bool PTUnrooted::isCanonicalName(const string& taxon, TaxonLevel level) {
	return StringUtils::startsWith(taxon, taxonLevel2prefix(level));
}

inline string PTUnrooted::PTUNode::getTaxon(double maxDist) const {
	return annoDist <= maxDist ? anno : anno + ";Other";
}

inline int8_t PTUnrooted::inferState(const Vector4d& loglik) {
	int8_t state = 0;
	loglik.maxCoeff(&state);
	return state;
}

inline Vector4d PTUnrooted::inferWeight(const Vector4d& loglik) {
	Vector4d p = (loglik.array() - loglik.maxCoeff()).exp(); /* scale before exponent */
	return p / p.sum();
}

inline Matrix4d PTUnrooted::initLeafMat() {
	Matrix4d leafMat = Matrix4d::Constant(infV);
	leafMat.diagonal().setConstant(0);
	return leafMat;
}

inline double PTUnrooted::estimateBranchLength(const Matrix4Xd& U, const Matrix4Xd& V,
		int start, int end, const string& method) {
	if(method == "unweighted")
		return estimateBranchLengthUnweighted(U, V, start, end);
	else if(method == "weighted")
		return estimateBranchLengthWeighted(U, V, start, end);
	else
		throw std::invalid_argument("Unknown branch length estimating method '" + method + "'");
}

inline ostream& PTUnrooted::PTPlacement::write(ostream& out) const {
	out << getId() << "\t" << ratio << "\t"
			<< getTaxonId() << "\t" << getTaxonName() << "\t"
			<< getAnnoDist() << "\t" << loglik << "\t"
			<< qPlace << "\t" << qTaxon;
	return out;
}

inline ostream& operator<<(ostream& out, const PTUnrooted::PTPlacement& place) {
	return place.write(out);
}

inline bool operator<(const PTUnrooted::PTLoc& lhs, const PTUnrooted::PTLoc& rhs) {
	return lhs.dist < rhs.dist;
}

inline bool compareByLoglik(const PTUnrooted::PTPlacement& lhs, const PTUnrooted::PTPlacement& rhs) {
	return lhs.loglik < rhs.loglik;
}

inline bool compareByQPlace(const PTUnrooted::PTPlacement& lhs, const PTUnrooted::PTPlacement& rhs) {
	return lhs.qPlace < rhs.qPlace;
}

inline bool compareByQTaxon(const PTUnrooted::PTPlacement& lhs, const PTUnrooted::PTPlacement& rhs) {
	return lhs.qTaxon < rhs.qTaxon;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREEUNROOTED_H_ */
