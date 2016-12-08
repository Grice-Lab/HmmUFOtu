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
#include <cmath>

namespace EGriceLab {
using namespace std;
using namespace EGriceLab;

const double PhyloTreeUnrooted::MAX_COST_EXP = -DBL_MIN_EXP / 2; /* use half of the DBL_MIN_EXP to avoid numeric-underflow */

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
			node2cost[*u][*v].setConstant(inf);
}

void PhyloTreeUnrooted::initInCost() {
	node2cost[root][root->parent] = Matrix4Xd::Constant(4, csLen, inf); /* initiate the dummy root -> NULL cost */
	for(vector<PTUNodePtr>::iterator u = id2node.begin(); u != id2node.end(); ++u)
		for(vector<PTUNodePtr>::iterator v = (*u)->neighbors.begin(); v != (*u)->neighbors.end(); ++v)
			node2cost[*u][*v] = Matrix4Xd::Constant(4, csLen, inf);
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
		for(vector<PTUNodePtr>::iterator child = node->neighbors.begin(); child != node->neighbors.end(); ++child) { /* check each child */
			if(!isChild(*child, node)) /* ignore non-child neighbor */
				continue;

			/* check whether this incoming cost from child has been evaluated */
			Vector4d inCost;
			if(isEvaluated(*child, node, j))
				inCost = node2cost[*child][node].col(j); /* use cached incoming cost */
			else {
				// cerr << "P:" << P << endl;
				inCost = evaluate(*child, j, model); /* evaluate this child recursively */
				if( (inCost.array() != inCost.array()).any() ) {
					cerr << "Warning: NaN values generated in cost from node id " << (*child)->id << " name " << (*child)->name
							<< " to node id " << node->id << " name " << node->name <<
							" with cost " << inCost.transpose() << endl;
				}
				node2cost[*child][node].col(j) = inCost; /* cache this cost */
			}

			/* convolute the inCost with Pr */
			double scale = 0; /* potential scaling factor to avoid numeric overflow */
			if((inCost.array() != inf).any() && (inCost.array() > MAX_COST_EXP).all()) { /* need re-scale only the min cost efficient is too large */
				scale = inCost.minCoeff() - MAX_COST_EXP;
				inCost.array() -= scale;
//				cerr << "scaling at " << scale << " to avoid numeric underflow" << endl;
			}
			const Matrix4d& P = model.Pr(node2length[*child][node]);

			for(Vector4d::Index i = 0; i < inCost.rows(); ++i)
				costVec(i) += -::log(P.row(i).dot((-inCost).array().exp().matrix())) + scale;

		}
//		cerr << "Internal node " << node->id << " evaluated at site " << j << " cost: "  << costVec.transpose() << endl;
	}

	if(node->isRoot()) /* the inCost message of root also needed to be stored */
		node2cost[node][node->parent].col(j) = costVec;
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

double PTUnrooted::treeCost(int j, const DNASubModel& model) {
	Vector4d cost = !isEvaluated(root, root->parent, j) ? /* tree site not evaluated */
		evaluate(root, j, model) /* evaluate if neccessary */ : node2cost[root][NULL].col(j) /* use cached value */;

	double scale = 0; /* scale factor to avoid potential numeric underflow */
	if((cost.array() != inf).any() && (cost.array() > MAX_COST_EXP).all()) { /* need re-scale only the min cost efficient is too large */
		scale = cost.minCoeff() - MAX_COST_EXP;
		cost.array() -= scale;
	}
	return -::log(model.getPi().dot((-cost).array().exp().matrix())) + scale;
}

} /* namespace EGriceLab */
