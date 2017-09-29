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
 * PhyloTreeUnrooted.cpp
 *
 *  Created on: Dec 1, 2016
 *      Author: zhengqi
 */

#include <sstream>
#include <stack>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <boost/algorithm/string.hpp> /* for boost string split */
#include <boost/lexical_cast.hpp>
#include "HmmUFOtuConst.h"
#include "ProgLog.h"
#include "StringUtils.h"
#include "SeqUtils.h"
#include "PhyloTreeUnrooted.h"
#include "DNASubModelFactory.h"

namespace EGriceLab {
using namespace std;
using namespace EGriceLab;
using Eigen::Map;
using Eigen::Matrix4Xd;

const double PhyloTreeUnrooted::MIN_LOGLIK_EXP = DBL_MIN_EXP / 2; /* use half of the DBL_MIN_EXP to avoid numeric-underflow */
const double PhyloTreeUnrooted::INVALID_LOGLIK = 1;
const double PhyloTreeUnrooted::LOGLIK_REL_EPS = 1e-6;
const double PhyloTreeUnrooted::BRANCH_EPS = 1e-5;

const string PhyloTreeUnrooted::KINDOM_PREFIX = "k__";
const string PhyloTreeUnrooted::PHYLUM_PREFIX = "p__";
const string PhyloTreeUnrooted::CLASS_PREFIX = "c__";
const string PhyloTreeUnrooted::ORDER_PREFIX = "o__";
const string PhyloTreeUnrooted::FAMILY_PREFIX = "f__";
const string PhyloTreeUnrooted::GENUS_PREFIX = "g__";
const string PhyloTreeUnrooted::SPECIES_PREFIX = "s__";
const Matrix4d PhyloTreeUnrooted::leafMat = initLeafMat();

static const char* TAXON_SEP = ";: "; /* valid taxon name separator */

bool PTUnrooted::isTip(const PTUNodePtr& node) {
	if(node->isLeaf())
		return false;
	for(vector<PTUNodePtr>::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child)
		if(isChild(*child, node) && !(*child)->isLeaf())
			return false;
	return true;
}

istream& PhyloTreeUnrooted::PhyloTreeUnrootedNode::load(istream& in) {
	/* read length */
	size_t nName, nAnno;
	in.read((char*) &nName, sizeof(string::size_type));
	in.read((char*) &nAnno, sizeof(string::size_type));

	/* read basic info */
	in.read((char*) &id, sizeof(long));
	StringUtils::loadString(name, in, nName);

	/* read seq */
	seq.load(in);

	/* read annotation */
	StringUtils::loadString(anno, in, nAnno);
	in.read((char*) &annoDist, sizeof(double));

	return in;
}

ostream& PhyloTreeUnrooted::PhyloTreeUnrootedNode::save(ostream& out) const {
	/* write length */
	size_t nName = name.length();
	size_t nAnno = anno.length();
	out.write((const char*) &nName, sizeof(size_t));
	out.write((const char*) &nAnno, sizeof(size_t));

	/* write basic info */
	out.write((const char*) &id, sizeof(long));
	StringUtils::saveString(name, out);

	/* write seq */
	seq.save(out);

	/* write annotation */
	StringUtils::saveString(anno, out);
	out.write((const char*) &annoDist, sizeof(double));

	return out;
}


PhyloTreeUnrooted::PhyloTreeUnrooted(const NewickTree& ntree) : csLen(0) {
	/* construct PTUNode by DFS of the NewickTree */
	boost::unordered_set<const NT*> visited;
	stack<const NT*> S;
	long id = 0; /* id start from 0 */
	unordered_map<const NT*, PTUNodePtr> nTree2PTree;

	S.push(&ntree);
	while(!S.empty()) {
		const NT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* construct this PTUNode */
			PTUNodePtr u = boost::make_shared<PTUNode>(id++, v->name);

			id2node.push_back(u);
			nTree2PTree[v] = u;

			/* add check each child of v */
			for(vector<NT>::const_iterator child = v->children.begin(); child != v->children.end(); ++child)
				S.push(&*child);
		}
	}

	/* explore the nTree again to establish the parent/child relationship */
	visited.clear();
	S.push(&ntree);
	while(!S.empty()) {
		const NT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* get corresponding PTUNode */
			const PTUNodePtr& u = nTree2PTree[v];
			if(root == NULL) // root node of the Newick tree encountered
				root = u;

			/* add check each child of u */
			for(vector<NT>::const_iterator Nchild = v->children.begin(); Nchild != v->children.end(); ++Nchild) {
				const PTUNodePtr& Pchild = nTree2PTree[&*Nchild];
				/* add this new edge */
				addEdge(u, Pchild);
				/* set parent */
				Pchild->parent = u;
				/* update branch length */
				setBranchLength(u, Pchild, Nchild->length);
				S.push(&*Nchild);
			}
		}
	}
	assert(root != NULL);
}

unsigned PhyloTreeUnrooted::loadMSA(const MSA& msa) {
	unsigned n0 = msaId2node.size(); /* original number of loaded nodes */
	if(msa.getAbc()->getAlias() != "DNA") {
		cerr << "PhyloTreeUnrooted can only read in MSA in DNA alphabet" << endl;
		return -1;
	}
	const unsigned numSeq = msa.getNumSeq();
	csLen = msa.getCSLen();

	/* check uniqueness of seq names in msa */
	unordered_map<string, unsigned> name2msaId;
	for(unsigned i = 0; i < numSeq; ++i) {
		string name = msa.seqNameAt(i);
		if(name2msaId.find(name) != name2msaId.end()) {
			cerr << "Non-unique seq name " << name << " found in your MSA data " << msa.getName() << endl;
			return -1;
		}
		else
			name2msaId[name] = i;
	}

	/* assign seq to each nodes of the tree, ignore nodes cannot be found (unnamed, etc) */
	for(vector<PTUNodePtr>::iterator node = id2node.begin(); node != id2node.end(); ++node) {
		assert(node - id2node.begin() == (*node)->id);

		unordered_map<string, unsigned>::const_iterator result = name2msaId.find((*node)->name);
		if(result == name2msaId.end()) /* this name cannot be found in the msa */
			continue;
		(*node)->seq = msa.dsAt(result->second);
		msaId2node[result->second] = *node;
		node2msaId[*node] = result->second;
	}
	assert(msaId2node.size() == node2msaId.size());
	return msaId2node.size() - n0;
}

istream& PTUnrooted::loadAnnotation(istream& in) {
	string line, name, anno;
	unordered_map<string, string> name2anno;
	while(getline(in, line)) {
		istringstream lineIn(line);
		std::getline(lineIn, name, ANNO_FIELD_SEP);
		std::getline(lineIn, anno, ANNO_FIELD_SEP);
		name2anno[name] = anno;
	}

	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		unordered_map<string, string>::const_iterator result = name2anno.find((*node)->name);
		if(result != name2anno.end())
			(*node)->name = result->second;
	}

	return in;
}

PhyloTreeUnrooted::PTUNodePtr PhyloTreeUnrooted::setRoot(const PTUNodePtr& newRoot) {
	if(newRoot == NULL || newRoot == root) /* no need to set */
		return root;

	newRoot->parent.reset(); // root has no parent
//	node2loglik[newRoot][NULL] = Matrix4Xd::Constant(4, csLen, inf); // new cache for dummy branch
	/* DFS of this tree starting from newRoot */
	boost::unordered_set<PTUNodePtr> visited;
	stack<PTUNodePtr> S;

	S.push(newRoot);
	while(!S.empty()) {
		PTUNodePtr u = S.top();
		S.pop();
		if(visited.find(u) == visited.end()) { /* not visited before */
			visited.insert(u);

			/* check each neighbor of v */
			for(vector<PTUNodePtr>::iterator v = u->neighbors.begin(); v != u->neighbors.end(); ++v) {
				if(visited.find(*v) == visited.end() /* not parent/ancestor of v */
					&& !isChild(*v, u)) /* update this child's parent */
					(*v)->parent = u;
				S.push(*v);
			}
		}
	}
	PTUNodePtr oldRoot = root;
	root = newRoot;
	return oldRoot;
}

void PhyloTreeUnrooted::calcNodeHeight() {
	for(vector<PTUNodePtr>::const_iterator leaf = id2node.begin(); leaf != id2node.end(); ++leaf) {
		if(!(*leaf)->isLeaf())
			continue;
		/* trace back this lineage */
		double h = 0;
		for(PTUnrooted::PTUNodePtr node = *leaf; node != NULL; node = node->parent) {
			if(node2height.find(node) == node2height.end() || h < node2height[node]) /* first time or shorter */
				node2height[node] = h;
			if(!node->isRoot())
				h += getBranchLength(node, node->getParent());
		}
	}
}

void PhyloTreeUnrooted::resetBranchLoglik() {
	for(vector<PTUNodePtr>::iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			node2branch[*u][*v].loglik.setConstant(INVALID_LOGLIK);
}

void PhyloTreeUnrooted::initBranchLoglik() {
	for(vector<PTUNodePtr>::iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v) /* u->neighbors */
			node2branch[*u][*v].loglik = Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK);
}

Vector4d PhyloTreeUnrooted::loglik(const PTUNodePtr& node, int j, double r) {
	Vector4d loglikVec = Vector4d::Zero();
	for(vector<PTUNodePtr>::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) {
		if(isChild(*child, node)) /* a child neighbor */
			loglikVec += dot_product_scaled(model->Pr(getBranchLength(node, *child) * r), loglik(*child, j)); /* evaluate child recursively */
	}
	if(node->isLeaf() && !node->seq.empty()) /* is a leaf root with assigned seq */
		loglikVec += getLeafLoglik(node->seq, j);
	return loglikVec;
}

Vector4d PhyloTreeUnrooted::loglik(const PTUNodePtr& node, int j) {
	if(isEvaluated(node, node->parent, j))
		return getBranchLoglik(node, node->parent, j);

	if(dG == NULL) {
		const Vector4d& loglikVec = loglik(node, j, 1); // using fixed rate
		/* cache this conditional loglik */
		setBranchLoglik(node, node->parent, j, loglikVec);
		return loglikVec;
	}
	else { /* use Gamma model */
		Matrix4Xd loglikMat(4, dG->getK());
		for(int k = 0; k < dG->getK(); ++k)
			loglikMat.col(k) = loglik(node, j, dG->rate(k));
		const Vector4d& loglikVec = row_mean_exp_scaled(loglikMat); // use average of DiscreteGammaModel rate
		/* cache this conditional loglik, parent could be NULL */
		setBranchLoglik(node, node->parent, j, loglikVec);
		return loglikVec;
	}
}

Matrix4Xd PTUnrooted::loglik(const PTUNodePtr& node) {
	if(isEvaluated(node, node->parent)) /* already evaluated */
		return getBranchLoglik(node, node->parent);

	Matrix4Xd loglikMat(4, csLen);
	for(int j = 0; j < csLen; ++j)
		loglikMat.col(j) = loglik(node, j);
	setBranchLoglik(node, node->parent, loglikMat);
	return loglikMat;
}


void PhyloTreeUnrooted::evaluate(const PTUNodePtr& node, int j) {
	if(isEvaluated(node, node->parent, j)) /* already evaluated */
		return;
	for(vector<PTUNodePtr>::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) { /* check each child */
		if(isChild(*child, node)) /* a child neighbor */
			loglik(*child, j); /* evaluate child recursively */
	}
}

ostream& PTUnrooted::exportTreeNewick(ostream& out, const PTUNodePtr& node) const {
	/* recursive write children */
	bool first = true;
	if(node->isRoot() || node->isInternal()) {
		out << '(';
		for(std::vector<PTUNodePtr>::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) {
			if(!isChild(*child, node)) /* not a child */
				continue;
			out << (first ? "" : ",");
			exportTreeNewick(out, *child);
			first = false;
		}
		out << ')';
	}
	/* write self */
	out << node->id;
	double length = getBranchLength(node, node->parent);
	if(length > 0)
		out << ':' << length;

	return out;
}

ostream& PTUnrooted::exportTreeNewick(ostream& out, const PTUNodePtr& node, const unordered_set<PTUNodePtr>& subset) const {
	/* test whether any child is flagged */
	bool flag = false;
	for(std::vector<PTUNodePtr>::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) {
		if(isChild(*child, node) && subset.find(*child) != subset.end()) {
			flag = true;
			break;
		}
	}
	if(flag) {
		bool first = true;
		if(node->isRoot() || node->isInternal()) {
			out << '(';
		}
		for(std::vector<PTUNodePtr>::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) {
			if(!isChild(*child, node)) /* not a child */
				continue;
			out << (first ? "" : ",");
			exportTreeNewick(out, *child, subset);
			first = false;
		}
		out << ')';
	}
	/* write self */
	out << node->id;
	double length = getBranchLength(node, node->parent);
	if(length > 0)
		out << ':' << length;

	return out;
}


vector<Matrix4d> PTUnrooted::getModelTraningSetGoldman() const {
	debugLog << "Training data using Gojobori method" << endl;
	vector<Matrix4d> data; // store observed base transition counts
	/* check every node of this tree */
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		if((*node)->isTip() && (*node)->neighbors.size() > 2) { // tip with >=2 children
//			cerr << "Find a candidate node id: " << (*node)->id << endl;
//			cerr << "First child: " << (*node)->firstChild()->name << endl;
//			cerr << "Last child: " << (*node)->lastChild()->name << endl;
			const DigitalSeq& seq1 = (*node)->firstChild()->seq;
			const DigitalSeq& seq2 = (*node)->lastChild()->seq;
			if(SeqUtils::pDist(seq1, seq1) <= DNASubModel::MAX_PDIST)
				data.push_back(DNASubModel::calcTransFreq2Seq(seq1, seq2));
		}
	}
	return data;
}

vector<Matrix4d> PTUnrooted::getModelTraningSetGojobori() const {
	vector<Matrix4d> data; // store observed base transition counts
	/* check every node of this tree */
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		const vector<PTUNodePtr> children = (*node)->getChildren();
		if(children.size() == 2 &&
				(children[0]->isTip() || children[1]->isTip()) ) { /* one child is a tip node */
			PTUNodePtr tipChild = children[0];
			PTUNodePtr outerChild = children[1];
			if(!tipChild->isTip())
				tipChild.swap(outerChild);

			const DigitalSeq& seq0 = PTUnrooted::randomLeaf(outerChild)->seq;
			const DigitalSeq& seq1 = tipChild->firstChild()->seq;
			const DigitalSeq& seq2 = tipChild->lastChild()->seq;
			if(SeqUtils::pDist(seq0, seq1) <= DNASubModel::MAX_PDIST &&
					SeqUtils::pDist(seq0, seq2) <= DNASubModel::MAX_PDIST)
								data.push_back(DNASubModel::calcTransFreq3Seq(seq0, seq1, seq2));
		}
	}
	debugLog << "Gojobori data prepared" << endl;
	return data;
}

Vector4d PTUnrooted::getModelFreqEst() const {
	Vector4d freq = Vector4d::Zero();
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		if((*node)->isLeaf())
			freq += DNASubModel::calcBaseFreq((*node)->seq);
	return freq;
}

istream& PTUnrooted::load(istream& in) {
	/* Read program info */
	string pname, pver;
	readProgName(in, pname);
	if(pname != progName) {
		errorLog << "Not a PTUnrooted object file" << endl;
		in.setstate(ios_base::failbit);
		return in;
	}
	readProgVersion(in, pver);
	if(cmpVersion(progVersion, pver) < 0) {
		errorLog << "You are trying using an older version " << (progName + progVersion) <<
				" to read a newer PTUnrooted data file that was build by " << (pname + pver) << endl;
		in.setstate(ios_base::failbit);
		return in;
	}

	/* init leaf matrix that does not depend on anything */
	initLeafMat();

	/* read global information */
	size_t nNodes;
	in.read((char*) &nNodes, sizeof(size_t));
	in.read((char*) &csLen, sizeof(int));

	/* read each node */
	for(size_t i = 0; i < nNodes; ++i) {
		PTUNodePtr node(new PTUNode); /* construct a new node */
		node->load(in);
		id2node.push_back(node);
	}

	/* read all edges */
	size_t nEdges;
	in.read((char*) &nEdges, sizeof(size_t));
	for(size_t i = 0; i < nEdges; ++i)
		loadEdge(in);

	/* load root */
	loadRoot(in);

	/* load node height */
	loadNodeHeight(in);

	/* load root loglik */
//	loadRootLoglik(in);

	/* read index */
	loadMSAIndex(in);

	/* load models */
	loadModel(in);
	loadDGModel(in);

	return in;
}

ostream& PTUnrooted::save(ostream& out) const {
	/* save program info */
	writeProgName(out, progName);
	writeProgVersion(out, progVersion);

	/* write global information */
	size_t nNodes = numNodes();
	out.write((const char*) &nNodes, sizeof(size_t));
	out.write((const char*) &csLen, sizeof(int));

	/* write each node */
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		(*node)->save(out);
	/* write all edges */
	size_t nEdges = numEdges();
	out.write((const char*) &nEdges, sizeof(size_t));
	for(vector<PTUNodePtr>::const_iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::const_iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			saveEdge(out, *u, *v);

	/* save root */
	saveRoot(out);

	/* save node height */
	saveNodeHeight(out);

	/* write index */
	saveMSAIndex(out);

	/* save models */
	saveModel(out);
	saveDGModel(out);

	return out;
}

ostream& PTUnrooted::saveMSAIndex(ostream& out) const {
	unsigned N = msaId2node.size();
	out.write((const char*) &N, sizeof(unsigned));
	for(map<unsigned, PTUNodePtr>::const_iterator it = msaId2node.begin(); it != msaId2node.end(); ++it) {
		out.write((const char*) &(it->first), sizeof(unsigned));
		out.write((const char*) &(it->second->id), sizeof(long));
	}

	return out;
}

istream& PTUnrooted::loadMSAIndex(istream& in) {
	unsigned N = 0;
	unsigned msaId;
	long id;
	in.read((char*) &N, sizeof(unsigned));
	for(unsigned i = 0; i < N; ++i) {
		in.read((char*) &msaId, sizeof(unsigned));
		in.read((char*) &id, sizeof(long));
		msaId2node[msaId] = id2node.at(id); /* build forward index */
		node2msaId[id2node.at(id)] = msaId; /* build reverse index */
	}

	return in;
}

ostream& PTUnrooted::saveEdge(ostream& out, const PTUNodePtr& node1, const PTUNodePtr& node2) const {
	out.write((const char*) &(node1->id), sizeof(long));
	out.write((const char*) &(node2->id), sizeof(long));
	bool flag = isParent(node1, node2);
	out.write((const char*) &flag, sizeof(bool));
	getBranch(node1, node2).save(out); /* save branch data */

	return out;
}

istream& PTUnrooted::loadEdge(istream& in) {
	long id1, id2;
	bool isParent;
	in.read((char*) &id1, sizeof(long));
	in.read((char*) &id2, sizeof(long));
	in.read((char*) &isParent, sizeof(bool));

	const PTUNodePtr& node1 = id2node[id1];
	const PTUNodePtr& node2 = id2node[id2];
	node1->neighbors.push_back(node2);
	if(isParent)
		node2->parent = node1;
	/* construct a new empty branch and load */
	node2branch[node1][node2].load(in);

	return in;
}

ostream& PTUnrooted::saveNodeHeight(ostream& out) const {
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		out.write((const char*) &((*node)->id), sizeof(long));
		out.write((const char*) &(node2height.at(*node)), sizeof(double));
	}

	return out;
}

istream& PTUnrooted::loadNodeHeight(istream& in) {
	size_t N = numNodes();
	long id;
	double h;
	for(size_t i = 0; i < N; ++i) {
		in.read((char*) &id, sizeof(long));
		in.read((char*) &h, sizeof(double));
		node2height[id2node[id]] = h;
	}

	return in;
}

istream& PTUnrooted::loadRoot(istream& in) {
	long rootId;
	double* buf = new double[4 * csLen];
	Map<Matrix4Xd> rootMap(buf, 4, csLen);
	/* set current root */
	in.read((char*) &rootId, sizeof(long));
	root = id2node[rootId];
	/* load current root loglik */
	in.read((char*) buf, 4 * csLen * sizeof(double));
	setBranchLoglik(root, NULL, rootMap);
	delete[] buf;

	return in;
}

ostream& PTUnrooted::saveRoot(ostream& out) const {
	/* save current root id */
	out.write((const char*) &(root->id), sizeof(long));
	double* buf = new double[4 * csLen];
	Map<Matrix4Xd> inLoglikMap(buf, 4, csLen);
	inLoglikMap = getBranchLoglik(root, NULL);
	out.write((const char*) buf, inLoglikMap.size() * sizeof(double));
	delete[] buf;

	return out;
}

istream& PTUnrooted::loadModel(istream& in) {
	string type, line;
	in >> type;
	in.ignore(); /* ignore the next '\n' character */
	/* create the model with a newly created object */
	model.reset(DNASubModelFactory::createModel(type));
	/* read model */
	in >> *model;
	return in;
}

ostream& PTUnrooted::saveModel(ostream& out) const {
	out << model->modelType() << endl;
	out << *model;
	return out;
}

istream& PTUnrooted::loadDGModel(istream& in) {
	bool modelSet;
	in.read((char*) &modelSet, sizeof(bool));
	if(modelSet) {
		dG.reset(new DiscreteGammaModel()); /* construct a new model and assign to dG */
		dG->load(in);
	}
	return in;
}

ostream& PTUnrooted::saveDGModel(ostream& out) const {
	bool modelSet = dG != NULL;
	out.write((const char*) &modelSet, sizeof(bool));
	if(modelSet)
		dG->save(out);
	return out;
}

PTUnrooted PTUnrooted::copySubTree(const PTUNodePtr& u, const PTUNodePtr& v) const {
	assert(isParent(v, u));

	PTUnrooted tree; /* construct an empty tree */
	long id = 0;
	tree.csLen = csLen; /* copy csLen */
	tree.model = model; /* copy the DNA model */
	tree.dG = dG; /* copy DiscreteGammaModel */

	/* construct new copy of nodes, but the old sequences are ignored */
	PTUNodePtr v2(new PTUNode(id++, v->name, v->anno, v->annoDist));
	PTUNodePtr u2(new PTUNode(id++, u->name, u->anno, u->annoDist));
	u2->parent = v2;

	/* add nodes */
	tree.id2node.push_back(v2);
	tree.id2node.push_back(u2);
	/* add edge */
	tree.addEdge(u2, v2);

	/* copy branch length and loglik */
	tree.setBranch(u2, v2, getBranch(u, v));
	tree.setBranch(v2, u2, getBranch(v, u));

	tree.setRoot(v2);
	return tree;
}

double PTUnrooted::optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v,
		int start, int end, double maxL) {
	assert(isParent(v, u));

	double w0 = getBranchLength(u, v);

	double q0 = ::exp(-w0);
	double p0 = 1 - q0;

	double p = p0;
	double q = q0;

	const Vector4d& pi = model->getPi();

	const Matrix4Xd& U = getBranchLoglik(u, v);
	const Matrix4Xd& V = getBranchLoglik(v, u);
	/* Felsenstein's iterative optimizing algorithm */
	for(int iter = 0; iter < MAX_ITER && p >= 0 && p <= 1; ++iter) {
		p = 0;
		int N = 0;
		for(int j = start; j <= end; ++j) {
			double logA = dot_product_scaled(pi, U.col(j) + V.col(j));
			double logB = dot_product_scaled(pi, U.col(j)) + dot_product_scaled(pi, V.col(j));
			if(::isnan(logA) || ::isnan(logB))
				continue;
			double scale = std::max(logA, logB);
			logA -= scale;
			logB -= scale;
			p += ::exp(logB) * p0 / (::exp(logA) * q0 + ::exp(logB) * p0);
			N++;
		}
		p /= N;
		q = 1 - p;

//		cerr << "N: " << N << " p: " << p << " q: " << q << endl;
		if(::fabs(::log(q) - ::log(q0)) < BRANCH_EPS)
			break;
		// update p0 and q0
		p0 = p;
		q0 = q;
	}

	double w = -::log(q); // final estimation
	if(w > maxL)
		w = maxL;
	setBranchLength(u, v, w);
//	cerr << "w0: " << w0 << " w: " << w << endl;

	return w;
}

double PTUnrooted::optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, const PTUNodePtr& r, const PTUNodePtr& n,
		int start, int end) {
	assert(root == r && isParent(r, u) && isParent(r, v) && isParent(r, n));

	double wur0 = getBranchLength(u, r);
	double wvr0 = getBranchLength(v, r);
	double wnr0 = getBranchLength(n, r);
	double w0 = wur0 + wvr0;

	double wur = wur0;
	double wvr = wvr0;
	double wnr = wnr0;

	/* every outgoing loglik(r,u), loglik(r,v) and loglik(r,n) depends on the other two incoming loglik */
	for(int iter = 0; iter < MAX_ITER && 0 <= wur && wur <= w0; ++iter) {
		/* evaluate loglik(r, n) and update wnr */
		setRoot(n);
		resetLoglik(r, n);
		evaluate(n, start, end);
		wnr = optimizeBranchLength(r, n, start, end, 1); /* do not use branch length > 1 */
		/* update loglik(r,u) and wur */
		setRoot(u);
		resetLoglik(r, u);
		evaluate(u, start, end);
		wur = optimizeBranchLength(r, u, start, end, w0);
		/* update wvr and loglik(r, v) */
		wvr = w0 - wur;
		setRoot(v);
		setBranchLength(r, v, wvr);
		resetLoglik(r, v);
		evaluate(v, start, end);

		setRoot(r);

		if(::abs(wur - wur0) < BRANCH_EPS && ::abs(wnr - wnr0) < BRANCH_EPS)
			break;

		wur0 = wur;
		wvr0 = wvr;
		wnr0 = wnr;
	}
//	cerr << "Estimated ratio: " << wur / w0 << endl;

	return wur / w0;
}

double PTUnrooted::estimateSeq(const DigitalSeq& seq, const PTUNodePtr& u, const PTUNodePtr& v,
		int start, int end, double ratio, double& wnr) const {
	assert(seq.length() == csLen);
	assert(isParent(v, u));

	double w0 = getBranchLength(u, v);
	const Matrix4Xd& U = getBranchLoglik(u, v);
	const Matrix4Xd& V = getBranchLoglik(v, u);
	const Matrix4Xd& N = getLeafLoglik(seq, start, end);

	double wur = w0 * ratio;
	double wvr = w0 * (1 - ratio);
	const Matrix4Xd& UPr = dot_product_scaled(model->Pr(wur), U, start, end); /* U*P(wur) */
	const Matrix4Xd& VPr = dot_product_scaled(model->Pr(wvr), V, start, end); /* V*P(wvr) */
	wnr = estimateBranchLength(UPr + VPr /* R */, N, start, end);

	return treeLoglik(model->getPi(), UPr + VPr + dot_product_scaled(model->Pr(wnr), N), /* N*P(wnr) */
				start, end);
}

double PTUnrooted::placeSeq(const DigitalSeq& seq, const PTUNodePtr& u, const PTUNodePtr& v,
		int start, int end, double ratio0, double wnr0) {
//	cerr << "Placing seq " << seq.getName() << " at " << u->getId() << "->" << v->getId() <<
//			" start: " << start << " end: " << end << " ratio0: " << ratio0 << " wnr0: " << wnr0 << endl;
	assert(seq.length() == csLen); /* make sure this is an aligned seq */
	assert(isParent(v, u));
	assert(0 <= ratio0 && ratio0 <= 1);

	/* break the connection of u and v */
	double w0 = getBranchLength(u, v);
	removeEdge(u, v);
	/* create a new interior root */
	PTUNodePtr r(new PTUNode(numNodes(), ""));
	/* create a new leaf with given seq */
	PTUNodePtr n(new PTUNode(numNodes() + 1, "", seq));
	n->parent = r;
	u->parent = r;
	v->parent = r;
	setRoot(r);
	/* add new nodes */
	id2node.push_back(r);
	id2node.push_back(n);
	/* place r at the ratio0 = wur0 / w0 */
	addEdge(u, r);
	addEdge(v, r);
	setBranch(u, r, getBranch(u, v));
	setBranch(v, r, getBranch(v, u));
	setBranchLength(u, r, w0 * ratio0);
	setBranchLength(v, r, w0 * (1 - ratio0));
	setBranchLoglik(r, u, Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK));
	setBranchLoglik(r, v, Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK));
	/* place r with initial branch length */
	addEdge(n, r);
	setBranchLength(n, r, wnr0);
	setBranchLoglik(r, n, Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK));
	setBranchLoglik(n, r, Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK));
	/* evaluate new incoming messages */
	evaluate(r, start, end); /* n->r evaluated */

	/* joint optimization */
	optimizeBranchLength(u, v, r, n, start, end);
	initRootLoglik();
	for(int j = start; j <= end; ++j) /* calculate root loglik */
		loglik(r, j);

	return treeLoglik(start, end);
}

bool PhyloTreeUnrooted::isFullCanonicalName(const string& taxon) {
	vector<string> fields;
	boost::split(fields, taxon, boost::is_any_of(TAXON_SEP), boost::token_compress_on);
	for(vector<string>::size_type level = 0; level < fields.size(); ++level)
		if(!isCanonicalName(fields[level], static_cast<TaxonLevel> (level)))
			return false;
	return true;
}

bool PhyloTreeUnrooted::isPartialCanonicalName(const string& taxon) {
	vector<string> fields;
	boost::split(fields, taxon, boost::is_any_of(TAXON_SEP), boost::token_compress_on);
	for(vector<string>::const_iterator name = fields.begin(); name != fields.end(); ++name)
		if(!isCanonicalName(*name))
			return false;
	return true;
}

string PhyloTreeUnrooted::formatTaxonName(const string& taxon) {
	if(taxon.empty())
		return taxon;

	vector<string> formatedTaxon;
	vector<string> fields;
	boost::split(fields, taxon, boost::is_any_of(TAXON_SEP), boost::token_compress_on);
	for(vector<string>::const_iterator name = fields.begin(); name != fields.end(); ++name)
		if(isCanonicalName(*name))
			formatedTaxon.push_back(*name);

	return boost::join(formatedTaxon, ";");
}

void PhyloTreeUnrooted::annotate() {
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		annotate(*node);
}

void PhyloTreeUnrooted::annotate(const PTUNodePtr& node) {
	vector<string> annoPath;
	PTUNodePtr p(node); /* pointer to current node */
	while(!isFullCanonicalName(p->name) && !p->isRoot()) { /* a non-full canonical named node */
		node->annoDist += getBranchLength(p, p->parent);
		if(isPartialCanonicalName(p->name))
			annoPath.push_back(p->name);
		p = p->parent;
	}
	if(isFullCanonicalName(p->name))
		annoPath.push_back(p->name); /* push last name */
	std::reverse(annoPath.begin(), annoPath.end()); /* reverse the annoPath */
	node->anno = !annoPath.empty() ? boost::join(annoPath, ";") : "Other";
}

size_t PhyloTreeUnrooted::estimateNumMutations(int j) const {
	size_t N = 0;
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		if(!(*node)->isRoot() && inferState((*node), j) != inferState((*node)->parent, j)) {
			N++;
		}
	}
	return N;
}

double PTUnrooted::estimateBranchLength(const Matrix4Xd& U, const Matrix4Xd& V, int start, int end) {
	assert(U.cols() == V.cols());
	assert(0 <= start && start <= end && end < U.cols());

	double d = 0;
	for(int j = start; j <= end; ++j) {
		const Vector4d& logU = U.col(j);
		const Vector4d& logV = V.col(j);
		int8_t b1 = inferState(logU);
		int8_t b2 = inferState(logV);
		if(b1 != b2) {
//			d++; /* unweighted count */
			d += inferWeight(logU)(b1) * inferWeight(logV)(b2); /* weighted count */
		}
	}
	return d / (end - start + 1);
}

double PTUnrooted::estimateBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, const DigitalSeq& seq,
		int start, int end) const {
	assert(isParent(v, u));
	assert(seq.length() == csLen);
	assert(0 <= start && start <= end && end < csLen);

	double w = getBranchLength(u, v);
	double d = 0;
	const Matrix4Xd& U = getBranchLoglik(u, v);
	const Matrix4Xd& V = getBranchLoglik(v, u);
	const Matrix4Xd& N = getLeafLoglik(seq, start, end);
	const Matrix4Xd& R = dot_product_scaled(model->Pr(w / 2), U, start, end) /* u->v convoluting to w/2 */
			           + dot_product_scaled(model->Pr(w / 2), V, start, end);

	for(int j = start; j <= end; ++j) {
		const Vector4d& logR = R.col(j);
		const Vector4d& logN = N.col(j);
		int8_t b1 = inferState(logR);
		int8_t b2 = inferState(logN);
		if(b1 != b2)
			d += inferWeight(logR)(b1) * inferWeight(logN)(b2); /* weighted count */
	}
	return d / (end - start + 1);
}

ostream& PTUnrooted::PTUBranch::save(ostream& out) const {
	out.write((const char*) &length, sizeof(double));
	size_t N = loglik.size();
	out.write((const char*) &N, sizeof(size_t));

	double *buf = new double[N];
	Map<Matrix4Xd> loglikMap(buf, 4, loglik.cols());
	loglikMap = loglik; /* copy data */
	out.write((const char*) buf, sizeof(double) * N);
	delete[] buf;

	return out;
}

istream& PTUnrooted::PTUBranch::load(istream& in) {
	in.read((char*) &length, sizeof(double));
	size_t N;
	in.read((char*) &N, sizeof(size_t));
	if(loglik.size() != N)
		loglik.resize(4, N / 4);

	double *buf = new double[N];
	in.read((char*) buf, sizeof(double) * N);
//	Map<Matrix4Xd> loglikMap(buf, 4, N / 4);
//	loglik = loglikMap;
	loglik = Map<Matrix4Xd>(buf, 4, N / 4); /* copy data */
	delete[] buf;

	return in;
}

void PTUnrooted::inferSeq(const PTUNodePtr& node) {
	if(node->seq.length() == csLen) /* already inferred */
		return;
	node->seq.setAbc(AlphabetFactory::nuclAbc); /* always use DNA alphabet */
	node->seq.resize(csLen);
	const Matrix4Xd& logMat = loglik(node);
	for(int j = 0; j < csLen; ++j)
		node->seq[j] = inferState(logMat.col(j));
}

DigitalSeq PTUnrooted::inferPostCS(const PTUNodePtr& node, const Matrix4Xd& count, double alpha) const {
	assert(count.cols() == csLen);
	/* construct the Dirichlet Prior */
	const Matrix4Xd& loglikMat = loglik(node);
	Matrix4Xd pri(4, csLen);
	for(int j = 0; j < csLen; ++j)
		pri.col(j) = inferWeight(loglikMat.col(j));
	Matrix4Xd postP = alpha * pri + count;
	postP.array().rowwise() /= postP.colwise().sum().array(); /* normalize postP by cols */
	/* infer consensus */
	DigitalSeq seq(AlphabetFactory::nuclAbc, boost::lexical_cast<string> (node->getId()));
	for(int j = 0; j < csLen; ++j)
		seq.push_back(inferState(postP.col(j)));
	return seq;
}

DigitalSeq PTUnrooted::inferPostCS(const PTUNodePtr& node, const Matrix4Xd& count, const RowVectorXd& gap, double alpha) const {
	assert(count.cols() == csLen || gap.cols() == csLen);
	/* construct the Dirichlet Prior */
	const Matrix4Xd& loglikMat = loglik(node);
	Matrix4Xd pri(4, csLen);
	for(int j = 0; j < csLen; ++j)
		pri.col(j) = inferWeight(loglikMat.col(j));
	Matrix4Xd postP = alpha * pri + count;
	postP.array().rowwise() /= postP.colwise().sum().array(); /* normalize postP by cols */
	/* infer consensus */
	DigitalSeq seq(AlphabetFactory::nuclAbc, boost::lexical_cast<string> (node->getId()));
	for(int j = 0; j < csLen; ++j)
		seq.push_back(count.col(j).sum() >= gap(j) ? inferState(postP.col(j)) : DegenAlphabet::GAP_BASE);
	return seq;
}

ostream& PTUnrooted::exportTree(ostream& out, const PTUNodePtr& node,
		boost::unordered_set<PTUNodePtr> subset,
		string format) const {
	StringUtils::toLower(format);
	/* expand subset to all their ancestors */
	for(boost::unordered_set<PTUNodePtr>::const_iterator it = subset.begin(); it != subset.end(); ++it)
		for(PTUNodePtr node = *it; node != NULL; node = node->parent)
			subset.insert(node);

	if(format == "newick")
		return exportTreeNewick(out, node, subset) << ";";
	else {
		errorLog << "Cannot write PTUnrooted tree, unknown tree format '" << format << "'" << std::endl;
		out.setstate(std::ios_base::failbit);
		return out;
	}
}

} /* namespace EGriceLab */


