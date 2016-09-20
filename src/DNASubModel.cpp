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

/** Estimate the substitution parameters using Goldman algorithm */
/*Matrix4d DNASubModel::estimateSubRateGoldman(const PhyloTree* tree) {
	assert(tree->isRooted());
	Matrix4d freq = Matrix4d::Zero();
	const set<const PhyloTree::PhyloTreeNode*>& leafSet = tree->dfsLeaves(); // get all nodes
	vector<const PhyloTree::PhyloTreeNode*> leaves(leafSet.begin(), leafSet.end()); // transfer into a list
	 pair-wise count of mutations
	for(vector<const PhyloTree::PhyloTreeNode*>::size_type i = 0; i < leaves.size() - 1; ++i)
		for(vector<const PhyloTree::PhyloTreeNode*>::size_type j = i + 1; j < leaves.size(); ++j)
			freq += updateParams2Seq(leaves[i], leaves[j]);
	return freq;
}*/

/** Estimate the substitution parameters using Gojobori algorithm */
/*Matrix4d DNASubModel::estimateSubRateGojobori(const PhyloTree* tree) {
	assert(tree->isRooted());
	Matrix4d freq = Matrix4d::Zero();
	const set<const PhyloTree::PhyloTreeNode*>& leafSet = tree->dfsLeaves(); // get all leaves
	for(set<const PhyloTree::PhyloTreeNode*>::const_iterator it = leafSet.begin(); it != leafSet.end(); ++it) {
		for(const PhyloTree::PhyloTreeNode* node = *it; !node->isRoot(); node = node->parent) {
			if(node == *it)
				continue;  ignore itself
			const set<const PhyloTree::PhyloTreeNode*>& siblingLeafSet = tree->dfsLeaves(tree->getSibling(node));
			vector<const PhyloTree::PhyloTreeNode*> siblingLeaves(siblingLeafSet.begin(), siblingLeafSet.end());
			assert(siblingLeaves.size() >= 2);  make sure there exists triples
			for(vector<const PhyloTree::PhyloTreeNode*>::size_type i = 0; i < siblingLeaves.size() - 1; ++i)
				for(vector<const PhyloTree::PhyloTreeNode*>::size_type j = i + 1; j < siblingLeaves.size(); ++j)
					freq += updateParams3Seq(*it, siblingLeaves[i], siblingLeaves[j]);
		}
	}
	return freq;
}*/

Matrix4d DNASubModel::calcTransFreq2Seq(const DigitalSeq& seq1, const DigitalSeq& seq2) {
	assert(abc == seq1.getAbc() && abc == seq2.getAbc());
	assert(seq1.length() == seq2.length());
	Matrix4d freq = Matrix4d::Zero();

	const DigitalSeq::size_type L = seq1.length();
	for(DigitalSeq::size_type i = 0; i < L; ++i)
		if(seq1.isSymbol(i) && seq2.isSymbol(i)) // both not a gap
			freq(seq1[i], seq2[i])++;
	return freq;
}

Matrix4d DNASubModel::calcTransFreq3Seq(const DigitalSeq& outer,
		const DigitalSeq& seq1, const DigitalSeq& seq2) {
	assert(abc == outer.getAbc() && abc == seq1.getAbc() && abc == seq2.getAbc());
	assert(outer.length() == seq1.length() && outer.length() == seq2.length());
	Matrix4d freq = Matrix4d::Zero();

	const DigitalSeq::size_type L = outer.length();
	for(DigitalSeq::size_type i = 0; i < L; ++i) {
		int8_t b0 = outer[i];
		int8_t b1 = seq1[i];
		int8_t b2 = seq2[i];
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

double DNASubModel::pDist(const DigitalSeq& seq1, const DigitalSeq& seq2) {
	assert(abc == seq1.getAbc() && abc == seq2.getAbc());
	assert(seq1.length() == seq2.length());
	int d = 0;
	int N = 0;
	for(DigitalSeq::size_type i = 0; i < seq1.length(); ++i) {
		int b1 = seq1[i];
		int b2 = seq2[i];
		if(b1 >= 0 && b2 >= 0) {
			N++;
			if(b1 != b2)
				d++;
		}
	}
	return static_cast<double>(d) / N;
}

Matrix4d DNASubModel::scale(Matrix4d Q, Vector4d pi, double mu = 1.0) {
	double beta = pi.dot(Q.diagonal());
	return Q / -beta * mu;
}


Matrix4d DNASubModel::logQfromP(Matrix4d P, bool reversible) {
	if(reversible)
		P = (P + P.transpose()) / 2.0;
	/* normalize P */
	for(Matrix4d::Index i = 0; i < P.rows(); ++i)
		P.row(i) /= P.row(i).sum();

	/* do matrix log by diagonalizable matrix decomposition */
	Vector4d lambda; /* eigen values of P */
	EigenSolver<Matrix4d> es(Q);
	if(es.info() != Success) {
		cerr << "Cannot perform EigenSolver on observed frequency data" << endl;
		abort();
	}
	lambda = es.eigenvalues();
	Matrix4d U = es.eigenvectors();
	Matrix4d U_ = U.inverse();
	Matrix4d X = P.asDiagonal();
	return U * X.array().log() * U_;
}

Matrix4d DNASubModel::constrainedQfromP(Matrix4d P, bool reversible) {
	if(reversible)
		P = (P + P.transpose()) / 2.0;
	Vector4d Z = P.rowwise().sum(); // normalizing constants
	Matrix4d Q = Matrix4d::Zero();
	/* set the elements */
	for(Matrix4d::Index i = 0; i < Q.rows(); ++i) {
		for(Matrix4d::Index j = 0; j < Q.cols(); ++j) {
			if(i != j) {
				Q(i, j) = P(i, j) / Z(i); /* non-diagonal */
				Q(i, i) -= Q(i, j);       /* diagonal */
			}
		}
	}
	return Q;
}

} /* namespace EGriceLab */
