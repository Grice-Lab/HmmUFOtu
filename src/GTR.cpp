/*
 * GTR.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: zhengqi
 */

#include <cassert>
#include <stack>
#include <set>
#include <cmath>
#include "GTR.h"
#include "StringUtils.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const string name = "GTR";


istream& GTR::read(istream& in) {
	string line, tag, value;
	while(in >> tag) {
		if(tag[0] == '#') { /* comment or header */
			std::getline(in, line); /* ignore the entire line */
			continue;
		}
		if(tag == "Type:") {
			in >> tag; // read in model type
			if(tag != modelType()) {
				cerr << "Unmatched Model Type!" << endl;
				cerr << "Trying to read in a " << tag << " model into a " << modelType() << " object" << endl;
				in.setstate(ios_base::badbit);
				return in;
			}
		}
		else if(tag == "pi:") {
			for(Vector4d::Index i = 0; i != pi.rows(); ++i)
				in >> pi(i);
		}
		else if(tag == "R:") {
			for(Vector4d::Index i = 0; i < R.rows(); ++i)
				for(Vector4d::Index j = 0; j < R.cols(); ++j)
					in >> R(i, j);
		}
		else {
			cerr << "Un-recognized line found in GTR Model file: tag: " << tag << endl << line << endl;
			in.setstate(ios_base::badbit);
			return in;
		}
	}
	setQfromParams();
	cerr << "Q set " << endl;
	return in;
}

ostream& GTR::write(ostream& out) const {
	out << "# DNA Substitution Model" << endl;
	out << "Type: " << modelType() << endl;
	out << "pi: " << pi.transpose().format(FULL_FORMAT) << endl;
	out << "R:" << endl << R.format(FULL_FORMAT) << endl;
//	out << "Q:" << endl << Q << endl;
	return out;
}

void GTR::trainParamsByDataset(const vector<Matrix4d>& P_vec, const Vector4d& f) {
	/* estimate pi using mean f */
	pi = f / f.sum();
//	assert(isValidFreq(pi));
	cerr << "pi estimated: " << pi.transpose() << endl;

	/* estimate Q from P_vec */
	Q.setZero();
	vector<Matrix4d>::size_type N = 0;
	for(vector<Matrix4d>::const_iterator it = P_vec.begin(); it != P_vec.end(); ++it) {
		const Matrix4d& Qv = constrainedQfromP(*it);
		if(isValidRate(Qv)) {
			N++;
			Q += scale(Qv);
		}
	}
	Q /= N;
	cerr << "Q estimated: " << endl << Q << endl;
	assert(isValidRate(Q));

	/* Decomposite R from Q, as Qij = Rij * pi(j) */
	for(Matrix4d::Index j = 0; j < R.cols(); ++j)
		R.col(j) = Q.col(j) / pi(j);
	R.diagonal().setZero(); /* set diagonal to zeros */

	/* average R to make it symmetric */
	R = (R + R.transpose()) / 2;
	/* reset Q */
	setQfromParams();
//	cerr << "Final data used in Q traning: " << N << endl;
//	cerr << "Final R estimated: " << endl << R << endl;
//	cerr << "Final Q estimated: " << endl << Q << endl;
}

void GTR::setQfromParams() {
	/* reset Q */
	assert(R.diagonal().sum() == 0);
	for(Matrix4d::Index j = 0; j < R.cols(); ++j)
		Q.col(j) = R.col(j) * pi(j);
	/* make Q valid */
	for(Matrix4d::Index i = 0; i < Q.rows(); ++i)
		Q(i, i) = - Q.row(i).sum();
	Q = scale(Q); /* re-scale Q */

	/* Eigen-decompsite Q
	 * Since Q is guaranteed to be a diagonizable matrix as Q = U * X * U-1
	 * Eigen values and vectors of Q are guaranteed to be real (not complex)
	 */
	EigenSolver<Matrix4d> es(Q);
	if(es.info() != Eigen::Success) {
		cerr << "Cannot perform EigenSolver on rate matrix Q:" << endl << Q << endl;
		abort();
	}
	lambda = es.eigenvalues().real();
	U = es.eigenvectors().real();
	U_1 = U.inverse();
}

void GTR::evaluate(PhyloTree& tree, int j) const {
	if(tree.isEvaluated(j)) // already evaluated
		return;
	/* calculate the cost of this node */
	if(tree.isLeaf()) { /* this is a leaf node */
		if(tree.seq[j] >= 0) { /* is not a gap */
			for(Matrix4Xd::Index s = 0; s < tree.cost.rows(); ++s)
				tree.cost(s, j) = s == tree.seq[j] ? 0 : inf;
		}
		else { /* a gap is treated as missing data */
			tree.cost.col(j) = - pi.array().log();
		}
//		if((tree.cost.col(j).array() != inf).count() != 1 && !(tree.cost.col(j).array() != inf).all())
//			cerr << "Strange Leaf " << tree.name << " cost j: " << tree.cost.col(j).transpose() << endl;
//		cerr << "Leaf " << tree.name << " evaluated. cost: \n" << tree.cost.col(j) << endl;
	}
	else { /* this is an internal node, all bases are treated as missing data */
		tree.cost.col(j).setZero(); /* initiate all cost to 0 */
		for(vector<PT>::iterator o = tree.children.begin(); o != tree.children.end(); ++o) { /* check each child */
			const Matrix4d& P = Pr(o->length);
//			cerr << "P:" << P << endl;
			evaluate(*o, j); /* evaluate this child recursively */
			for(Matrix4Xd::Index i = 0; i < tree.cost.rows(); ++i)
				tree.cost(i, j) += -::log(P.row(i).dot((-o->cost.col(j)).array().exp().matrix())) + o->scale(j);
		}
		/* rescale the node if neccessary */
		double maxC = tree.cost.col(j).maxCoeff();
		if(maxC >= PhyloTree::MIN_EXPONENT)
			tree.reScale(maxC - PhyloTree::MIN_EXPONENT, j);
//
//		cerr << "Internal node evaluated " << tree.name << " evaluated. cost: \n" << tree.cost.col(j) << endl;
	}

}

} /* namespace EGriceLab */
