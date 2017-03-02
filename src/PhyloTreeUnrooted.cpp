/*
 * PhyloTreeUnrooted.cpp
 *
 *  Created on: Dec 1, 2016
 *      Author: zhengqi
 */

#include <sstream>
#include <stack>
#include <boost/unordered_set.hpp>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <algorithm>
#include "HmmUFOtuConst.h"
#include "ProgLog.h"
#include "StringUtils.h"
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
const double PhyloTreeUnrooted::BRANCH_EPS = 1e-6;

const string PhyloTreeUnrooted::KINDOM_PREFIX = "k__";
const string PhyloTreeUnrooted::PHYLUM_PREFIX = "p__";
const string PhyloTreeUnrooted::CLASS_PREFIX = "c__";
const string PhyloTreeUnrooted::ORDER_PREFIX = "o__";
const string PhyloTreeUnrooted::FAMILY_PREFIX = "f__";
const string PhyloTreeUnrooted::GENUS_PREFIX = "g__";
const string PhyloTreeUnrooted::SPECIES_PREFIX = "s__";

bool PTUnrooted::isTip(const PTUNodePtr& node) {
	if(node->isLeaf())
		return false;
	for(BranchMap::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child)
		if(isChild(child->first, node) && !child->first->isLeaf())
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


PhyloTreeUnrooted::PhyloTreeUnrooted(const NewickTree& nTree) : csLen(0) {
	/* construct PTUNode by DFS of the NewickTree */
	boost::unordered_set<const NT*> visited;
	stack<const NT*> S;
	long id = 0; /* id start from 0 */
	boost::unordered_map<const NT*, PTUNodePtr> n2pTree;

	S.push(&nTree);
	while(!S.empty()) {
		const NT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* construct this PTUNode */
			PTUNodePtr u(new PTUNode(id++, v->name));

			id2node.push_back(u);
			n2pTree[v] = u;

			/* add each child of v */
			for(vector<NT>::const_iterator nChild = v->children.begin(); nChild != v->children.end(); ++nChild)
				S.push(&*nChild);
		}
	}

	/* explore the nTree again to establish the parent/child relationship */
	visited.clear();
	S.push(&nTree);
	while(!S.empty()) {
		const NT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* get corresponding PTUNode */
			const PTUNodePtr& u = n2pTree[v];
			if(root == NULL)
				root = u;

			/* add check each child of u */
			for(vector<NT>::const_iterator nChild = v->children.begin(); nChild != v->children.end(); ++nChild) {
				const PTUNodePtr& pChild = n2pTree[&*nChild];
				/* update parent add pChild as a neighbor of parent */
				u->neighbors[pChild] = PTUBranch(nChild->length);
				/* update child by setting its parent */
				pChild->neighbors[u] = PTUBranch(nChild->length);
				pChild->parent = u;

				S.push(&*nChild);
			}
		}
	}

}

size_t PhyloTreeUnrooted::loadMSA(const MSA& msa) {
	const DegenAlphabet* abc = msa.getAbc();
	if(abc->getAlias() != "DNA") {
		cerr << "PhyloTreeUnrooted can only read in MSA in DNA alphabet" << endl;
		return -1;
	}
	const unsigned numSeq = msa.getNumSeq();
	csLen = msa.getCSLen();

	/* check uniqueness of seq names in msa */
	map<string, unsigned> nameIdx;
	for(unsigned i = 0; i != numSeq; ++i) {
		string name = msa.seqNameAt(i);
		if(nameIdx.find(name) != nameIdx.end()) {
			cerr << "Non-unique seq name " << name << " found in your MSA data " << msa.getName() << endl;
			return -1;
		}
		else {
			nameIdx[name] = i;
		}
	}
	size_t assigned = 0;
	/* assign seq to each nodes of the tree, ignore nodes cannot be found (unnamed, etc) */
	for(vector<PTUNodePtr>::iterator node = id2node.begin(); node != id2node.end(); ++node) {
		assert(node - id2node.begin() == (*node)->id);

		map<string, unsigned>::const_iterator result = nameIdx.find((*node)->name);
		if(result == nameIdx.end()) /* this name cannot be found in the msa */
			continue;
		(*node)->seq = msa.dsAt(result->second);
		assigned++;
	}
	return assigned;
}

istream& PTUnrooted::loadAnnotation(istream& in) {
	string line, name, anno;
	boost::unordered_map<string, string> name2anno;
	while(getline(in, line)) {
		istringstream lineIn(line);
		std::getline(lineIn, name, ANNO_FIELD_SEP);
		std::getline(lineIn, anno, ANNO_FIELD_SEP);
		name2anno[name] = anno;
	}

	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		boost::unordered_map<string, string>::const_iterator result = name2anno.find((*node)->name);
		if(result != name2anno.end())
			(*node)->name = result->second;
	}

	return in;
}

PhyloTreeUnrooted::PTUNodePtr PhyloTreeUnrooted::setRoot(const PTUNodePtr& newRoot) {
	if(newRoot == NULL || newRoot == root) /* no need to set */
		return root;

	newRoot->parent = NULL; // root has no parent
//	node2loglik[newRoot][NULL] = Matrix4Xd::Constant(4, csLen, inf); // new cache for dummy branch
	/* DFS of this tree starting from newRoot */
	boost::unordered_set<PTUNodePtr> visited;
	stack<PTUNodePtr> S;

	S.push(newRoot);
	while(!S.empty()) {
		PTUNodePtr v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);

			/* check each neighbor of v */
			for(BranchMap::iterator neighbor = v->neighbors.begin(); neighbor != v->neighbors.end(); ++neighbor) {
				if(visited.find(neighbor->first) == visited.end() /* not visited, not parent/ancestor of v */
					&& !isChild(neighbor->first, v)) /* not updated yet */
					neighbor->first->parent = v;
				S.push(neighbor->first);
			}
		}
	}
	PTUNodePtr oldRoot = root;
	root = newRoot;
	return oldRoot;
}


void PhyloTreeUnrooted::resetLoglik() {
	for(vector<PTUNodePtr>::iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(BranchMap::iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			v->second.loglik.setConstant(INVALID_LOGLIK);
}

void PhyloTreeUnrooted::initInLoglik() {
	for(vector<PTUNodePtr>::iterator u = id2node.begin(); u != id2node.end(); ++u) {
		for(BranchMap::iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			v->second.loglik = Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK);
	}
}

void PhyloTreeUnrooted::initLeafLoglik() {
	leafLoglik.resize(4, 5);
	if(model == NULL) /* no model provided yet */
		leafLoglik.setConstant(INVALID_LOGLIK);
	else {
		/* initiate the non-gap (0..3) columns */
		leafLoglik.leftCols(4).setConstant(infV);
		leafLoglik.leftCols(4).diagonal().setConstant(0);
		leafLoglik.col(4) = model->getPi().array().log();
	}
}

Vector4d PhyloTreeUnrooted::loglik(const PTUNodePtr& node, int j, double r) {
	Vector4d loglikVec = Vector4d::Zero();
	for(BranchMap::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) {
		if(isChild(child->first, node)) /* a child neighbor */
			loglikVec += dot_product_scaled(model->Pr(child->second.length * r), loglik(child->first, j)); /* evaluate child recursively */
	}
	if(node->isLeaf() && !node->seq.empty()) /* is a leaf root with assigned seq */
		loglikVec += node->seq[j] >= 0 ? leafLoglik.col(node->seq[j]) /* a base observed */ : leafLoglik.col(4) /* a gap observed */;
	return loglikVec;
}

Vector4d PhyloTreeUnrooted::loglik(const PTUNodePtr& node, int j) {
	if(isEvaluated(node, node->parent, j))
		return getBranchLoglik(node, node->parent, j);

	Vector4d loglikVec;
	if(dG == NULL)
		loglikVec = loglik(node, j, 1); // using fixed rate
	else {
		Matrix4Xd loglikMat(4, dG->getK());
		for(int k = 0; k < dG->getK(); ++k)
			loglikMat.col(k) = loglik(node, j, dG->rate(k));
		loglikVec = row_mean_exp_scaled(loglikMat); // use average of DiscreteGammaModel rate
	}
	/* cache this conditional loglik for non-root node */
	if(!node->isRoot())
		setBranchLoglik(node, node->parent, j, loglikVec);
	return loglikVec;
}

Matrix4Xd PTUnrooted::loglik(const PTUNodePtr& node) {
	if(isEvaluated(node, node->parent)) /* not a root and evaluated */
		return getBranchLoglik(node, node->parent);

	Matrix4Xd loglikMat(4, csLen);
	for(int j = 0; j < csLen; ++j)
		loglikMat.col(j) = loglik(node, j);
	if(!node->isRoot())
		setBranchLoglik(node, node->parent, loglikMat);
	return loglikMat;
}

double PTUnrooted::treeLoglik(const PTUNodePtr& node, int start, int end) {
	double loglik = 0;
	for(int j = start; j <= end; ++j)
		loglik += treeLoglik(node, j);
	return loglik;
}

int PhyloTreeUnrooted::inferState(const PTUnrooted::PTUNodePtr& node, int j) {
	if(!node->seq.empty())
		return node->seq[j];
	Vector4d logV = loglik(node, j);
	int state;
	logV.maxCoeff(&state);
	return state;
}


void PhyloTreeUnrooted::evaluate(const PTUNodePtr& node, int j) {
	if(isEvaluated(node, node->parent, j)) /* already evaluated */
		return;
	for(BranchMap::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) { /* check each child */
		if(isChild(child->first, node)) /* a child neighbor */
			loglik(child->first, j); /* evaluate child recursively */
	}
}

ostream& PTUnrooted::writeTreeNewick(ostream& out, const PTUNodePtr& node) const {
	bool first = true;
	if(node->isRoot() || node->isInternal()) {
		out << '(';
		for(BranchMap::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) {
			if(!isChild(child->first, node)) /* not a child */
				continue;
			out << (first ? "" : ",");
			writeTreeNewick(out, child->first);
			first = false;
		}
		out << ')';
	}
	if(StringUtils::containsWhiteSpace(node->name) || StringUtils::containsAny(node->name, NewickTree::INVALID_CHARS)) // name contains INVALID CHARS
		out << "'" << node->name << "'";
	else
		out << node->name;
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
			if(DNASubModel::pDist(seq1, seq1) <= DNASubModel::MAX_PDIST)
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
				(children.front()->isTip() || children.back()->isTip()) ) { /* one child is a tip node */
//			cerr << "Find a candidate node id: " << (*node)->id << endl;
//			cerr << "Child 1 is tip " << children.front()->isTip() << endl;
//			cerr << "Child 2 is tip " << children.back()->isTip() << endl;
			PTUNodePtr tipChild = children.front();
			PTUNodePtr outerChild = children.back();
			if(!tipChild->isTip())
				tipChild.swap(outerChild);

			const DigitalSeq& seq0 = PTUnrooted::randomLeaf(outerChild)->seq;
			const DigitalSeq& seq1 = tipChild->firstChild()->seq;
			const DigitalSeq& seq2 = tipChild->lastChild()->seq;
			if(DNASubModel::pDist(seq0, seq1) <= DNASubModel::MAX_PDIST &&
					DNASubModel::pDist(seq0, seq2) <= DNASubModel::MAX_PDIST)
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
//	initInLoglik();
//	initLeafLoglik();
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

	/* read leaf loglik */
	loadLeafLoglik(in);

	/* load root */
	loadRoot(in);

	/* load root loglik */
//	loadRootLoglik(in);

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
		for(BranchMap::const_iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			saveEdge(out, *u, v->first, v->second);

	/* write leaf loglik */
	saveLeafLoglik(out);
//	debugLog << "leaf loglik saved" << endl;

	/* save root */
	saveRoot(out);
//	debugLog << "Root saved" << endl;

	/* save models */
	saveModel(out);
	saveDGModel(out);

	return out;
}

ostream& PTUnrooted::saveEdge(ostream& out, const PTUNodePtr& u, const PTUNodePtr& v, const PTUBranch& w) const {
	out.write((const char*) &(u->id), sizeof(long));
	out.write((const char*) &(v->id), sizeof(long));
	bool flag = isParent(u, v);
	out.write((const char*) &flag, sizeof(bool));
	w.save(out);
	return out;
}

istream& PTUnrooted::loadEdge(istream& in) {
	long id1, id2;
	bool isParent;
	in.read((char*) &id1, sizeof(long));
	in.read((char*) &id2, sizeof(long));
	in.read((char*) &isParent, sizeof(bool));

	PTUNodePtr u = id2node[id1];
	PTUNodePtr v = id2node[id2];
	if(isParent)
		v->parent = u;

	u->neighbors[v].load(in); /* add a new branch u->v, then read its data */

	return in;
}

istream& PTUnrooted::loadLeafLoglik(istream& in) {
	leafLoglik.resize(4, 5);
	double* buf = new double[leafLoglik.size()];
	in.read((char*) buf, leafLoglik.size() * sizeof(double));
	leafLoglik = Map<Matrix4Xd>(buf, leafLoglik.rows(), leafLoglik.cols()); /* copy via Map */
	delete[] buf;
	return in;
}

ostream& PTUnrooted::saveLeafLoglik(ostream& out) const {
	double* buf = new double[leafLoglik.size()];
	Map<Matrix4Xd> leafLoglikMap(buf, leafLoglik.rows(), leafLoglik.cols());
	leafLoglikMap = leafLoglik; /* copy via Map */
	out.write((const char*) buf, leafLoglik.size() * sizeof(double));
	delete[] buf;
	return out;
}

istream& PTUnrooted::loadRoot(istream& in) {
	/* set current root */
	long rootId;
	in.read((char*) &rootId, sizeof(long));
	root = id2node[rootId];

	return in;
}

ostream& PTUnrooted::saveRoot(ostream& out) const {
	/* save current root id */
	out.write((const char*) &(root->id), sizeof(long));
	return out;
}

istream& PTUnrooted::loadRootLoglik(istream& in) {
	long id;
	double* buf = new double[4 * csLen];
	Map<Matrix4Xd> inLoglikMap(buf, 4, csLen);

	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		in.read((char*) &id, sizeof(long));
		assert(id == (*node)->id);
		in.read((char*) buf, 4 * csLen * sizeof(double));
		/* assign root loglik */
		(*node)->neighbors[NULL].loglik = inLoglikMap;
	}
	delete[] buf;

	return in;
}

ostream& PTUnrooted::saveRootLoglik(ostream& out) const {
	double* buf = new double[4 * csLen];
	Map<Matrix4Xd> inLoglikMap(buf, 4, csLen);

	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		out.write((const char*) &((*node)->id), sizeof(long));
		inLoglikMap = getBranchLoglik(*node, NULL); /* copy node->NULL loglik */
		out.write((const char*) buf, inLoglikMap.size() * sizeof(double));
	}
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
	assert(isChild(u, v));

	PTUnrooted tree; /* construct an empty tree */
	long id = 0;
	tree.csLen = csLen; /* copy csLen */
	tree.model = model; /* copy the DNA model */
	tree.dG = dG; /* copy DiscreteGammaModel */
	tree.leafLoglik = leafLoglik; /* copy leaf loglik */

	/* construct new copies of nodes */
	PTUNodePtr v2(new PTUNode(id++, v->name, v->seq, v->anno, v->annoDist));
	PTUNodePtr u2(new PTUNode(id++, u->name, u->seq, u->anno, u->annoDist));
	u2->neighbors[v2] = u->neighbors[v]; // copy u->v branch length and loglik
	v2->neighbors[u2] = v->neighbors[u]; // copy v->u branch length and loglik
	u2->parent = v2;

	/* add nodes into subtree */
	tree.id2node.push_back(v2);
	tree.id2node.push_back(u2);

	tree.setRoot(v2);

	return tree;
}

double PTUnrooted::estimateBranchLength(const PTUNodePtr& u, const PTUNodePtr& v, int start, int end) {
	assert(isParent(v, u));
	const Vector4d& pi = model->getPi();

	const Matrix4Xd& loglikU = getBranchLoglik(u, v);
	const Matrix4Xd& loglikV = getBranchLoglik(v, u);

	double p = 0;
	for(int j = start; j <= end; ++j) {
		double logA = dot_product_scaled(pi, loglikU.col(j) + loglikV.col(j));
		double logB = dot_product_scaled(pi, loglikU.col(j)) + dot_product_scaled(pi, loglikV.col(j));
		double scale = std::max(logA, logB);
		logA -= scale;
		logB -= scale;
		p += ::exp(logB) / (::exp(logA) + ::exp(logB));
	}
	return p / (end - start + 1);
}

double PTUnrooted::optimizeBranchLength(const PTUNodePtr& u, const PTUNodePtr& v,
		int start, int end) {
	assert(isParent(v, u));

	double w0 = getBranchLength(u, v);

	double q0 = ::exp(-w0);
	double p0 = 1 - q0;

	double p = p0;
	double q = q0;

	const Vector4d& pi = model->getPi();

	const Matrix4Xd& loglikU = getBranchLoglik(u, v);
	const Matrix4Xd& loglikV = getBranchLoglik(v, u);
	/* Felsenstein's iterative optimizing algorithm */
	while(p >= 0 && p <= 1) {
		p = 0;
		for(int j = start; j <= end; ++j) {
			double logA = dot_product_scaled(pi, loglikU.col(j) + loglikV.col(j));
			double logB = dot_product_scaled(pi, loglikU.col(j)) + dot_product_scaled(pi, loglikV.col(j));
			double scale = std::max(logA, logB);
			logA -= scale;
			logB -= scale;
			p += ::exp(logB) * p0 / (::exp(logA) * q0 + ::exp(logB) * p0);
		}
		p /= (end - start + 1);
		q = 1 - p;

		if(::fabs(q - q0) < BRANCH_EPS)
			break;
		// update p0 and q0
		p0 = p;
		q0 = q;
	}

	double w = -::log(q0); // final estimation
	setBranchLength(u, v, w);
	return w;
}

PTUnrooted& PTUnrooted::placeSeq(const DigitalSeq& seq, const PTUNodePtr& u, const PTUNodePtr& v,
		int start, int end) {
	assert(seq.length() == csLen); /* make sure this is an aligned seq */
	assert(isParent(v, u));

	/* break the connection of u and v */
	u->neighbors.erase(v);
	v->neighbors.erase(u);

	/* create a new interior root */
	PTUNodePtr r(new PTUNode(numNodes(), v->name, csLen));
	/* create a new leaf with given seq */
	PTUNodePtr n(new PTUNode(numNodes() + 1, v->name, seq));
	id2node.push_back(r);
	id2node.push_back(n);

	/* place r at the middle-point of u and v */
	double w0 = getBranchLength(u, v);
	addEdge(u, r, u->neighbors[v], w0 / 2); /* u->r mimic original u->v with half length */
	addEdge(r, u, w0 / 2, Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK)); /* r->u not evaluted */
	addEdge(v, r, v->neighbors[u], w0 / 2); /* v->r mimic original v->u with half length */
	addEdge(r, v, w0 / 2, Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK)); /* r->v not evaluated */

	/* place new node n */
	addEdge(n, r, 0, Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK));
	addEdge(r, n, 0, Matrix4Xd::Constant(4, csLen, INVALID_LOGLIK));

	/* evaluate subtree */
	setRoot(n);
	evaluate(n); /* r->n loglik evaluated */
	setRoot(r);
	evaluate(r); /* n->r loglik evaluated */
	setBranchLength(n, r, estimateBranchLength(n, r, start, end)); /* n->r branch length estimated */

	double lnr = optimizeBranchLength(n, r, start, end);
//	double lur = optimizeBranchLength(u, r, start, end);
//	double lvr = node2length[r][v] = node2length[v][r] = v0 - vur;

	/* annotate the new nodes */
	r->anno = v->anno;
	r->annoDist = getBranchLength(r, v);
	n->anno = r->anno;
	n->annoDist = r->annoDist + lnr;
//	n->annoDist = node2length[r][n] / 2;

	return *this;
}

string& PhyloTreeUnrooted::formatTaxaName(string& taxa) {
	if(taxa.empty())
		return taxa;
	/* remove white spaces */
	taxa.erase(std::remove_if(taxa.begin(), taxa.end(), ::isspace), taxa.end());

	/* append tailing ; if necessary */
	if(taxa.back() != ';')
		taxa.push_back(';');

	/* remove empty taxa */
	StringUtils::removeAll(taxa, KINDOM_PREFIX + ';');
	StringUtils::removeAll(taxa, PHYLUM_PREFIX + ';');
	StringUtils::removeAll(taxa, CLASS_PREFIX + ';');
	StringUtils::removeAll(taxa, ORDER_PREFIX + ';');
	StringUtils::removeAll(taxa, FAMILY_PREFIX + ';');
	StringUtils::removeAll(taxa, GENUS_PREFIX + ';');
	StringUtils::removeAll(taxa, SPECIES_PREFIX + ';');

	/* remove tailing ; */
	taxa.erase(taxa.end() - 1);

	return taxa;
}

vector<PTUnrooted::PTUNodePtr> PTUnrooted::getLeafHits(const DigitalSeq& seq, double maxPDist,
		int start, int end) const {
	vector<PTUnrooted::PTUNodePtr> hits;
	for(vector<PTUnrooted::PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
		if((*node)->isLeaf() && DNASubModel::pDist((*node)->getSeq(), seq, start, end + 1) <= maxPDist)
			hits.push_back(*node);
	}
	return hits;
}

void PhyloTreeUnrooted::annotate() {
	for(vector<PTUNodePtr>::const_iterator nodeIt = id2node.begin(); nodeIt != id2node.end(); ++nodeIt) {
		/* find first named ancestor of this node, if any */
		PTUNodePtr node(*nodeIt);
		while(!node->isRoot() && !isCanonicalName(node->name)) {
			(*nodeIt)->annoDist += getBranchLength(node, node->parent);
			node = node->parent;
		}
		(*nodeIt)->anno = isCanonicalName(node->name) ?
				*nodeIt == node ? node->name : node->name + ";Other"
						: "Other";
	}
}

size_t PhyloTreeUnrooted::estimateNumMutations(int j) {
	size_t N = 0;
	for(vector<PTUNodePtr>::const_iterator nodeIt = id2node.begin(); nodeIt != id2node.end(); ++nodeIt) {
		if(!(*nodeIt)->isRoot() && inferState((*nodeIt), j) != inferState((*nodeIt)->parent, j)) {
			N++;
		}
	}
	return N;
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

	double *buf = new double[N];
	in.read((char*) buf, sizeof(double) * N);
//	Map<Matrix4Xd> loglikMap(buf, 4, N / 4);
//	loglik = loglikMap;
	loglik = Map<Matrix4Xd>(buf, 4, N / 4); /* copy data */
	delete[] buf;

	return in;
}

} /* namespace EGriceLab */


