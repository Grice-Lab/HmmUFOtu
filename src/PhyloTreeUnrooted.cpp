/*
 * PhyloTreeUnrooted.cpp
 *
 *  Created on: Dec 1, 2016
 *      Author: zhengqi
 */

#include <stack>
#include <boost/unordered_set.hpp>
#include <cfloat>
#include <cmath>
#include "HmmUFOtuConst.h"
#include "PhyloTreeUnrooted.h"
#include "DNASubModelFactory.h"

namespace EGriceLab {
using namespace std;
using namespace EGriceLab;
using Eigen::Map;
using Eigen::Matrix4Xd;

const double PhyloTreeUnrooted::MAX_COST_EXP = -DBL_MIN_EXP / 2; /* use half of the DBL_MIN_EXP to avoid numeric-underflow */
const double PhyloTreeUnrooted::INVALID_COST = -1;

istream& PhyloTreeUnrooted::PhyloTreeUnrootedNode::load(istream& in) {
	char* buf = NULL;
	string::size_type nName, nAnno;
	DigitalSeq::size_type nSeq;

	/* read basic info */
	in.read((char*) &id, sizeof(long));
	in.read((char*) &nName, sizeof(string::size_type));
	buf = new char[nName + 1];
	in.read(buf, nName + 1); /* read the null terminal */
	name.assign(buf, nName); // override the original value
	delete[] buf;

	in.read((char*) &nSeq, sizeof(DigitalSeq::size_type));
	buf = new char[nSeq + 1];
	in.read(buf, nSeq + 1); /* read the null terminal */
	seq.assign((const int8_t*) buf, nSeq);
	delete[] buf;

	in.read((char*) &nAnno, sizeof(string::size_type));
	buf = new char[nAnno + 1];
	in.read(buf, nAnno + 1); /* read the null terminal */
	anno.assign(buf, nAnno);
	delete[] buf;

	in.read((char*) &annoDist, sizeof(double));

	return in;
}

ostream& PhyloTreeUnrooted::PhyloTreeUnrootedNode::save(ostream& out) const {

	/* get aux length */
	string::size_type nName = name.length();
	DigitalSeq::size_type nSeq = seq.length();
	string::size_type nAnno = anno.length();
	/* write basic info */
	out.write((const char*) &id, sizeof(long));
	out.write((const char*) &nName, sizeof(string::size_type));
	out.write(name.c_str(), nName + 1);
	out.write((const char*) &nSeq, sizeof(DigitalSeq::size_type));
	out.write((const char*) seq.c_str(), nSeq + 1);
	out.write((const char*) &nAnno, sizeof(string::size_type));
	out.write(anno.c_str(), nAnno + 1);
	out.write((const char*) &annoDist, sizeof(double));
	return out;
}


PhyloTreeUnrooted::PhyloTreeUnrooted(const NewickTree& ntree) : csLen(0) {
	/* construct PTUNode by DFS of the NewickTree */
	boost::unordered_set<const NT*> visited;
	stack<const NT*> S;
	long id = 0; /* id start from 0 */
	boost::unordered_map<const NT*, PTUNodePtr> nTree2PTree;

	S.push(&ntree);
	while(!S.empty()) {
		const NT* v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			visited.insert(v);
			/* construct this PTUNode */
			PTUNodePtr u(new PTUNode(id++, v->name));

			id2node.push_back(u);
			nTree2PTree[v] = u;

			/* add check each child of v */
			for(vector<NT>::const_iterator childIt = v->children.begin(); childIt != v->children.end(); ++childIt)
				S.push(&*childIt);
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
			if(root == NULL)
				root = u;

			/* add check each child of u */
			for(vector<NT>::const_iterator childIt = v->children.begin(); childIt != v->children.end(); ++childIt) {
				const PTUNodePtr& child = nTree2PTree[&*childIt];
				/* update parent */
				u->neighbors.push_back(child);
				/* update child */
				child->neighbors.push_back(u);
				child->parent = u;
				/* update branch length */
				node2length[u][child] = node2length[child][u] = childIt->length;

				S.push(&*childIt);
			}
		}
	}

}

long PhyloTreeUnrooted::loadMSA(const MSA& msa) {
	const DegenAlphabet* abc = msa.getAbc();
	if(abc != SeqCommons::nuclAbc) {
		cerr << "PhyloTreeUnrooted can only read in MSA in dna alphabet" << endl;
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
	long assigned = 0;
	/* assign seq to each nodes of the tree, ignore nodes cannot be found (unnamed, etc) */
	for(vector<PTUNodePtr>::iterator nodeIt = id2node.begin(); nodeIt != id2node.end(); ++nodeIt) {
		assert(nodeIt - id2node.begin() == (*nodeIt)->id);

		map<string, unsigned>::const_iterator result = nameIdx.find((*nodeIt)->name);
		if(result == nameIdx.end()) /* this name cannot be found in the msa */
			continue;
		(*nodeIt)->seq = msa.dsAt(result->second);
		assigned++;
	}
	return assigned;
}

PhyloTreeUnrooted::PTUNodePtr PhyloTreeUnrooted::setRoot(const PTUNodePtr& newRoot) {
	if(newRoot == NULL || newRoot == root) /* no need to set */
		return root;

	newRoot->parent = NULL; // root has no parent
//	node2cost[newRoot][NULL] = Matrix4Xd::Constant(4, csLen, inf); // new cache for dummy branch
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
			for(vector<PTUNodePtr>::iterator neighborIt = v->neighbors.begin(); neighborIt != v->neighbors.end(); ++neighborIt) {
				if(visited.find(*neighborIt) == visited.end() /* not parent/ancestor of v */
					&& !isChild(*neighborIt, v)) /* update this child's parent */
					(*neighborIt)->parent = v;
				S.push(*neighborIt);
			}
		}
	}
	PTUNodePtr oldRoot = root;
	root = newRoot;
	return oldRoot;
}

void PhyloTreeUnrooted::resetCost() {
	for(vector<PTUNodePtr>::iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			node2cost[*u][*v].setConstant(INVALID_COST);
}

void PhyloTreeUnrooted::initInCost() {
	for(vector<PTUNodePtr>::iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v) /* u->children */
			node2cost[*u][*v] = Matrix4Xd::Constant(4, csLen, INVALID_COST);
}

void PhyloTreeUnrooted::initLeafCost() {
	leafCost.resize(4, 5);
	if(model == NULL) /* no model provided yet */
		leafCost.setConstant(INVALID_COST);
	else {
		/* initiate the non-gap (0..3) columns */
		leafCost.leftCols(4).setConstant(inf);
		leafCost.leftCols(4).diagonal().setConstant(0);
		leafCost.col(4) = - model->getPi().array().log();
	}
}

Vector4d PhyloTreeUnrooted::evaluate(const PTUNodePtr& u, const PTUNodePtr& v, int j) {
	if(isEvaluated(u, v, j))
		return node2cost[u][v].col(j);
	Vector4d costVec;

	if(u->isLeaf()) { /* this is a leaf node */
		costVec = u->seq[j] >= 0 ? leafCost.col(u->seq[j]) /* a base observed */ : leafCost.col(4) /* a gap observed */;
//		cerr << "Leaf node " << node->name << " evaluated at site " << j << " cost: "  << costVec.transpose() << endl;
	}
	else { /* this is an internal node, all bases are treated as missing data */
		costVec.setZero();
		for(vector<PTUNodePtr>::const_iterator child = u->neighbors.begin(); child != u->neighbors.end(); ++child) { /* check each child */
			if(!isChild(*child, u)) /* ignore non-child neighbor */
				continue;

			/* check whether this incoming cost from child has been evaluated */
			const Vector4d& inCost = evaluate(*child, u, j); /* evaluate recursively */
//				if( (inCost.array() != inCost.array()).any() ) {
//					errorLog << "Warning: NaN values generated in cost from node id " << (*child)->id << " name " << (*child)->name
//							 << " to node id " << node->id << " name " << node->name <<
//							" with cost " << inCost.transpose() << endl;
//				}
			/* convolute the inCost with Pr */
			costVec += dot_product_scaled(model->Pr(node2length[*child][u]), inCost);
		}
//		cerr << "Internal node " << node->id << " evaluated at site " << j << " cost: "  << costVec.transpose() << endl;
	}

	node2cost[u][v].col(j) = costVec;
	return costVec;
}

ostream& PTUnrooted::writeTreeNewick(ostream& out, const PTUNodePtr& node) const {
	bool first = true;
	if(node->isRoot() || node->isInternal()) {
		out << '(';
		for(std::vector<PTUNodePtr>::const_iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) {
			if(!isChild(*child, node)) /* not a child */
				continue;
			out << (first ? "" : ",");
			writeTreeNewick(out, *child);
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

double PTUnrooted::treeCost(int j) {
//	evaluate(root, j, model); /* evaluate first */
	Vector4d rootCost = Vector4d::Zero();
	for(vector<PTUNodePtr>::const_iterator child = root->neighbors.begin(); child != root->neighbors.end(); ++child)
		if(isChild(*child, root)) /* a real child */
			rootCost += dot_product_scaled(model->Pr(node2length[*child][root]), node2cost[*child][root].col(j));

	/* all incoming cost collected */
	if(root->isLeaf()) /* this is a leaf root, need add in the cost of leaf without convolution */
		rootCost += (root->seq[j] >= 0 ? leafCost.col(root->seq[j]) /* a base observed */ : leafCost.col(4)); /* a gap observed */

	/* final convolution */
	return dot_product_scaled(model->getPi(), rootCost);
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
//	initInCost();
//	initLeafCost();
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

	/* read edge costs */
	for(size_t i = 0; i < nEdges; ++i)
		loadEdgeCost(in);

	/* read leaf cost */
	loadLeafCost(in);

	/* load root */
	loadRoot(in);

	/* load model */
	loadModel(in);

	return in;
}

ostream& PTUnrooted::save(ostream& out) const {
	/* save program info */
	writeProgName(out, progName);
	writeProgVersion(out, progVersion);
//	debugLog << "program info saved" << endl;

	/* write global information */
	size_t nNodes = numNodes();
	out.write((const char*) &nNodes, sizeof(size_t));
	out.write((const char*) &csLen, sizeof(int));
//	debugLog << "global information saved" << endl;

	/* write each node */
	for(vector<PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node)
		(*node)->save(out);
//	debugLog << "all nodes saved" << endl;

	/* write all edges */
	size_t nEdges = numEdges();
	out.write((const char*) &nEdges, sizeof(size_t));
	for(vector<PTUNodePtr>::const_iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::const_iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			saveEdge(out, *u, *v);
//	debugLog << "all edge saved" << endl;

	/* write edge costs */
	for(vector<PTUNodePtr>::const_iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::const_iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			saveEdgeCost(out, *u, *v);
//	debugLog << "edge costs saved" << endl;

	/* write leaf cost */
	saveLeafCost(out);
//	debugLog << "leaf cost saved" << endl;

	/* save root */
	saveRoot(out);
//	debugLog << "Root saved" << endl;

	/* save model */
	saveModel(out);
//	debugLog << "Model saved" << endl;

	return out;
}

ostream& PTUnrooted::saveEdge(ostream& out, const PTUNodePtr& node1, const PTUNodePtr& node2) const {
	out.write((const char*) &(node1->id), sizeof(long));
	out.write((const char*) &(node2->id), sizeof(long));
	bool flag = isParent(node1, node2);
	out.write((const char*) &flag, sizeof(bool));
	double length = getBranchLength(node1, node2);
	out.write((const char*) &length, sizeof(double));

	return out;}

istream& PTUnrooted::loadEdge(istream& in) {
	long id1, id2;
	bool isParent;
	double length;
	in.read((char*) &id1, sizeof(long));
	in.read((char*) &id2, sizeof(long));
	in.read((char*) &isParent, sizeof(bool));
	in.read((char*) &length, sizeof(double));

	const PTUNodePtr& node1 = id2node[id1];
	const PTUNodePtr& node2 = id2node[id2];
	node1->neighbors.push_back(node2);
	if(isParent)
		node2->parent = node1;
	node2length[node1][node2] = length;

	return in;
}

istream& PTUnrooted::loadLeafCost(istream& in) {
	leafCost.resize(4, 5);
	double* buf = new double[leafCost.size()];
	in.read((char*) buf, leafCost.size() * sizeof(double));
	leafCost = Map<Matrix4Xd>(buf, leafCost.rows(), leafCost.cols()); /* copy via Map */
	delete[] buf;
	return in;
}

ostream& PTUnrooted::saveLeafCost(ostream& out) const {
	double* buf = new double[leafCost.size()];
	Map<Matrix4Xd> leafCostMap(buf, leafCost.rows(), leafCost.cols());
	leafCostMap = leafCost; /* copy via Map */
	out.write((const char*) buf, leafCost.size() * sizeof(double));
	delete[] buf;
	return out;
}

istream& PTUnrooted::loadEdgeCost(istream& in) {
	long id1, id2;
	in.read((char*) &id1, sizeof(long));
	in.read((char*) &id2, sizeof(long));

	double* buf = new double[4 * csLen];
	in.read((char*) buf, 4 * csLen * sizeof(double));
	Map<Matrix4Xd> inCostMap(buf, 4, csLen); /* a customized map */
	/* assign this cost */
	node2cost[id2node[id1]][id2node[id2]] = inCostMap;
	delete[] buf;
	return in;
}

ostream& PTUnrooted::saveEdgeCost(ostream& out, const PTUNodePtr& node1, const PTUNodePtr& node2) const {
	out.write((const char*) &(node1->id), sizeof(long));
	out.write((const char*) &(node2->id), sizeof(long));
	double* buf = new double[4 * csLen];
	Map<Matrix4Xd> inCostMap(buf, 4, csLen);
	inCostMap = getBranchCost(node1, node2); /* copy the matrix */
	out.write((const char*) buf, inCostMap.size() * sizeof(double));
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
	out << *model << endl;

	return out;
}

} /* namespace EGriceLab */


