/*
 * DNASubModel.cpp
 *
 *  Created on: Apr 1, 2016
 *      Author: zhengqi
 */

#include <vector>
#include <set>
#include <vector>
#include <stack>
#include "DNASubModel.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const double DNASubModel::MAX_PDIST = 0.15; /* maximum p-dist between training sequences */
const IOFormat DNASubModel::FULL_FORMAT(Eigen::FullPrecision);
const IOFormat DNASubModel::STD_FORMAT(Eigen::StreamPrecision);

Matrix4d DNASubModel::calcTransFreq2Seq(const DigitalSeq& seq1, const DigitalSeq& seq2) {
	assert(seq1.getAbc() == seq2.getAbc());
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
	assert(outer.getAbc() == seq1.getAbc() && outer.getAbc() == seq2.getAbc());
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

Vector4d DNASubModel::calcBaseFreq(const DigitalSeq& seq) {
	Vector4d f = Vector4d::Zero();
	for(DigitalSeq::const_iterator it = seq.begin(); it != seq.end(); ++it)
		if(*it >= 0)
			f(*it)++;
	return f;
}

double DNASubModel::pDist(const DigitalSeq& seq1, const DigitalSeq& seq2) {
	assert(seq1.getAbc() == seq2.getAbc());
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

Matrix4d DNASubModel::scale(Matrix4d Q, Vector4d pi, double mu) {
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
	Vector4cd lambda; /* eigen values of P */
	EigenSolver<Matrix4d> es(P);
	if(es.info() != Eigen::Success) {
		cerr << "Cannot perform EigenSolver on observed frequency data P:" << endl << P << endl;
		abort();
	}
	lambda = es.eigenvalues();
	Matrix4cd U = es.eigenvectors();
	Matrix4cd U_1 = U.inverse();
	Matrix4cd X = lambda.asDiagonal();
	return (U * X.array().log().matrix() * U_1).real();
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
//	cerr << "P: " << P << endl;
//	cerr << "Q: " << Q << endl;
	return Q;
}

void DNASubModel::trainParamsGoldman(const PhyloTree& tree) {
	vector<Matrix4d> P_vec; /* store observed base transition counts */
	Vector4d f; /* observed base counts */
	/* do DFS explore of the tree */
	stack<const PT*> S;
	set<const PT*> visited;
	S.push(&tree);
	while(!S.empty()) {
		const PT* node = S.top();
		S.pop();
		if(visited.find(node) == visited.end()) { /* node not visited */
			visited.insert(node);
			if(node->isTip() && node->children.size() >= 2) { /* use its first 2 children for training */
				const DigitalSeq& seq1 = node->children[0].seq;
				const DigitalSeq& seq2 = node->children[1].seq;
				if(DNASubModel::pDist(seq1, seq1) <= MAX_PDIST) {
					P_vec.push_back(DNASubModel::calcTransFreq2Seq(seq1, seq2));
					f += DNASubModel::calcBaseFreq(seq1);
					f += DNASubModel::calcBaseFreq(seq2);
				}
			}
			for(vector<PT>::const_iterator it = node->children.begin(); it != node->children.end(); ++it) { // its children
				const PT* p = &*it;
				S.push(p);
			}
		}
	}
	trainParamsByDataset(P_vec, f);
}

void DNASubModel::trainParamsGojobori(const PhyloTree& tree) {
	vector<Matrix4d> P_vec; /* store observed base transition counts */
	Vector4d f = Vector4d::Zero(); /* store observed base counts */
	/* do DFS explore of the tree */
	stack<const PT*> S;
	set<const PT*> visited;
	S.push(&tree);
	while(!S.empty()) {
		const PT* node = S.top();
		S.pop();
		if(visited.find(node) == visited.end()) { /* node not visited */
			visited.insert(node);
			if(node->children.size() == 2 &&
				(node->children[0].isTip() || node->children[1].isTip()) ) { /* one child is a tip node */
				const DigitalSeq& seq0 = node->children[0].isTip() ? node->children[1].firstLeaf()->seq : node->children[0].firstLeaf()->seq;
				const DigitalSeq& seq1 = node->children[0].isTip() ? node->children[0].children[0].seq : node->children[1].children[0].seq;
				const DigitalSeq& seq2 = node->children[0].isTip() ? node->children[0].children[1].seq : node->children[1].children[1].seq;

//				cerr << "node: " << node->name << endl;
//				cerr << "seq0.name : " << seq0.getName() << endl;
//				cerr << "seq1.name : " << seq1.getName() << endl;
//				cerr << "seq2.name : " << seq2.getName() << endl;

				if(DNASubModel::pDist(seq0, seq1) <= MAX_PDIST && DNASubModel::pDist(seq0, seq2) <= MAX_PDIST) {
					P_vec.push_back(DNASubModel::calcTransFreq3Seq(seq0, seq1, seq2));
					f += calcBaseFreq(seq0);
					f += calcBaseFreq(seq1);
					f += calcBaseFreq(seq2);
				}
			}
			for(vector<PT>::const_iterator it = node->children.begin(); it != node->children.end(); ++it) { // its children
				const PT* p = &*it;
				S.push(p);
			}
		}
	}
	trainParamsByDataset(P_vec, f);
}


} /* namespace EGriceLab */
