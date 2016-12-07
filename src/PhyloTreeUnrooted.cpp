/*
 * PhyloTreeUnrooted.cpp
 *
 *  Created on: Dec 1, 2016
 *      Author: zhengqi
 */

#include <stack>
#include <boost/unordered_set.hpp>
#include "PhyloTreeUnrooted.h"
#include <cfloat>

namespace EGriceLab {
using namespace std;
using namespace EGriceLab;

const double PhyloTreeUnrooted::MAX_COST_EXP = ::fabs(DBL_MIN_EXP) - ::fabs(FLT_MIN_EXP) > ::fabs(FLT_MIN_EXP) ? ::fabs(DBL_MIN_EXP) - ::fabs(FLT_MIN_EXP) : ::fabs(FLT_MIN_EXP);

PhyloTreeUnrooted::PhyloTreeUnrooted(const NewickTree& ntree) : L(0) {
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
	const unsigned csLen = msa.getCSLen();
	L = csLen;

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
	cerr << "Resetting root at " << newRoot->name << endl;

	newRoot->parent = NULL; // root has no parent
	/* DFS of this tree starting from newRoot */
	boost::unordered_set<PTUNodePtr> visited;
	stack<PTUNodePtr> S;

	S.push(newRoot);
	while(!S.empty()) {
		PTUNodePtr v = S.top();
		S.pop();
		if(visited.find(v) == visited.end()) { /* not visited before */
			cerr << "  setRoot is exploring " << v->name << endl;
			visited.insert(v);

			/* check each neighbor of v */
			for(vector<PTUNodePtr>::iterator neighborIt = v->neighbors.begin(); neighborIt != v->neighbors.end(); ++neighborIt) {
				cerr << "    setRoot found a neighbor of " << v->name << " as " << (*neighborIt)->name << endl;
				if(visited.find(*neighborIt) == visited.end()) { /* not parent/ancestor of v */
					/* update this child's parent */
					cerr << "    setting " << (*neighborIt)->name << " 's parent as " << v->name << endl;
					(*neighborIt)->parent = v;
				}
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
			node2cost[*u][*v].setConstant(inf);
}

void PhyloTreeUnrooted::initInCost() {
	for(vector<PTUNodePtr>::iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			node2cost[*u][*v] = Matrix4Xd::Constant(4, L, inf);
}

void PhyloTreeUnrooted::initLeafCost(const DNASubModel& model) {
	leafCost.resize(4, 5);
	/* initiate the non-gap (0..3) columns */
	leafCost.leftCols(4).setConstant(inf);
	leafCost.leftCols(4).diagonal().setConstant(0);
	leafCost.col(4) = - model.getPi().array().log();
}

Vector4d PhyloTreeUnrooted::evaluate(const PTUNodePtr& node, int j, const DNASubModel& model) {
	Vector4d costVec = Vector4d::Zero();

	if(node->isLeaf()) { /* this is a leaf node, always evaluated on the fly */
		return node->seq[j] >= 0 ? leafCost.col(node->seq[j]) /* a base observed */ : leafCost.col(4) /* a gap observed */;
//		cerr << "Leaf node " << node->name << " evaluated at site " << j << " cost: "  << costVec.transpose() << endl;
	}
	else { /* this is an internal node, all bases are treated as missing data */
		for(vector<PTUNodePtr>::iterator o = node->neighbors.begin(); o != node->neighbors.end(); ++o) { /* check each child */
			if(!isChild(*o, node)) /* ignore non-child neighbor */
				continue;

			/* check whether this incoming cost from child o has been evaluated */
			Vector4d inCost;
			if(isEvaluated(*o, node, j))
				inCost = node2cost[*o][node].col(j); /* use cached value */
			else {
				const Matrix4d& P = model.Pr(node2length[*o][node]);
				// cerr << "P:" << P << endl;
				Vector4d oCost = evaluate(*o, j, model); /* evaluate this child recursively */
				if((oCost.array() == inf).all()) {
					cerr << "strange ocost from node id " << (*o)->id << " name " << (*o)->name
							<< " to node id " << node->id << " name " << node->name <<
							" with cost " << oCost.transpose() << endl;
				}
				double scale = 0; /* potential scaling factor to avoid numeric overflow */
				if((oCost.array() > MAX_COST_EXP).all()) { /* need re-scale only the min cost efficient is too large */
					scale = oCost.minCoeff() - MAX_COST_EXP;
					oCost.array() -= scale;

					cerr << "scaling at " << scale << " to avoid numeric underflow" << endl;
				}
				for(Vector4d::Index i = 0; i < inCost.rows(); ++i)
					inCost(i) = -::log(P.row(i).dot((-oCost).array().exp().matrix())) + scale;

				node2cost[*o][node].col(j) = inCost; /* cache this cost */
			}
			costVec += inCost;
		}
//		cerr << "Internal node " << node->id << " evaluated at site " << j << " cost: "  << costVec.transpose() << endl;
	}

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

} /* namespace EGriceLab */
