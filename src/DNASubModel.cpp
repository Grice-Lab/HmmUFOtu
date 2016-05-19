/*
 * DNASubModel.cpp
 *
 *  Created on: Apr 1, 2016
 *      Author: zhengqi
 */

#include <vector>
#include <set>
#include <vector>
#include "DNASubModel.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

void DNASubModel::updateRate() {
	Q = R * pi;
	for(Matrix4d::Index i = 0; i < R.rows(); ++i) {
		double q = 0;
		for(Matrix4d::Index j = 0; j < R.cols(); ++j)
			if(i != j)
				q -= Q(i, j);
		Q(i, i) = q;
	}
}

void DNASubModel::scaleRate() {
	double beta = pi.dot(Q.diagonal());
	Q /= -beta;
}

/** Estimate the substitution parameters using Goldman algorithm */
Matrix4d DNASubModel::estimateSubRateGoldman(const PhyloTree* tree) {
	assert(tree->isRooted());
	Matrix4d freq = Matrix4d::Zero();
	const set<const PhyloTree::PhyloTreeNode*>& leafSet = tree->dfsLeaves(); // get all nodes
	vector<const PhyloTree::PhyloTreeNode*> leaves(leafSet.begin(), leafSet.end()); // transfer into a list
	/* pair-wise count of mutations */
	for(vector<const PhyloTree::PhyloTreeNode*>::size_type i = 0; i < leaves.size() - 1; ++i)
		for(vector<const PhyloTree::PhyloTreeNode*>::size_type j = i + 1; j < leaves.size(); ++j)
			freq += updateParams2Seq(leaves[i], leaves[j]);
	return freq;
}

/** Estimate the substitution parameters using Gojobori algorithm */
Matrix4d DNASubModel::estimateSubRateGojobori(const PhyloTree* tree) {
	assert(tree->isRooted());
	Matrix4d freq = Matrix4d::Zero();
	const set<const PhyloTree::PhyloTreeNode*>& leafSet = tree->dfsLeaves(); // get all leaves
	for(set<const PhyloTree::PhyloTreeNode*>::const_iterator it = leafSet.begin(); it != leafSet.end(); ++it) {
		for(const PhyloTree::PhyloTreeNode* node = *it; !node->isRoot(); node = node->parent) {
			if(node == *it)
				continue; /* ignore itself */
			const set<const PhyloTree::PhyloTreeNode*>& siblingLeafSet = tree->dfsLeaves(tree->getSibling(node));
			vector<const PhyloTree::PhyloTreeNode*> siblingLeaves(siblingLeafSet.begin(), siblingLeafSet.end());
			assert(siblingLeaves.size() >= 2); /* make sure there exists triples */
			for(vector<const PhyloTree::PhyloTreeNode*>::size_type i = 0; i < siblingLeaves.size() - 1; ++i)
				for(vector<const PhyloTree::PhyloTreeNode*>::size_type j = i + 1; j < siblingLeaves.size(); ++j)
					freq += updateParams3Seq(*it, siblingLeaves[i], siblingLeaves[j]);
		}
	}
	return freq;
}

Matrix4d DNASubModel::updateParams2Seq(const PhyloTree::PhyloTreeNode* seq1, const PhyloTree::PhyloTreeNode* seq2) {
	assert(seq1->seq.getAbc() == seq2->seq.getAbc() || *(seq1->seq.getAbc()) == *(seq2->seq.getAbc()));
	assert(seq1->seq.length() == seq2->seq.length());
	Matrix4d freq = Matrix4d::Zero();

	const DigitalSeq::size_type L = seq1->seq.length();
	for(DigitalSeq::size_type i = 0; i < L; ++i)
		if(seq1->seq[i] >= 0 && seq2->seq[i] >= 0) // both not a gap
			freq(seq1->seq[i], seq2->seq[i])++;
	return freq;
}

Matrix4d DNASubModel::updateParams3Seq(const PhyloTree::PhyloTreeNode* outer,
		const PhyloTree::PhyloTreeNode* seq1, const PhyloTree::PhyloTreeNode* seq2) {
	assert(outer->seq.getAbc() == seq1->seq.getAbc() && outer->seq.getAbc() == seq2->seq.getAbc() ||
			*(outer->seq.getAbc()) == *(seq1->seq.getAbc()) && *(outer->seq.getAbc()) == *(seq2->seq.getAbc()));
	assert(outer->seq.length() == seq1->seq.length() && outer->seq.length() == seq2->seq.length());
	Matrix4d freq = Matrix4d::Zero();

	const DigitalSeq::size_type L = outer->seq.length();
	for(DigitalSeq::size_type i = 0; i < L; ++i) {
		int8_t b0 = outer->seq[i];
		int8_t b1 = seq1->seq[i];
		int8_t b2 = seq2->seq[i];
		int8_t bc; // ancestor of b0, b1 and b2
		if(!(b0 >= 0 && b1 >= 0 &&  b2 >= 0)) // ignore any gaps
			continue;
		if(b0 == b1 && b0 == b2) // no change
			bc = b0;
		else if(b0 == b1 && b0 != b2) // change is between bc->b2
			bc = b0;
		else if(b0 == b2 && b0 != b1) // change is between bc->b1
			bc = b0;
		else if(b0 != b1 && b0 != b2 && b1 == b2) // change is between bc->b0
			bc = b1;
		else // all different, cannot guess
			continue;
		freq(bc, b0)++;
		freq(bc, b1)++;
		freq(bc, b2)++;
	}
	return freq;
}

} /* namespace EGriceLab */
