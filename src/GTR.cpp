/*
 * GTR.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: zhengqi
 */

#include <cassert>
#include <stack>
#include <set>
#include "GTR.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const string name = "GTR";

Matrix4d GTR::Pr(double t, double r) const {
	Matrix4d X = Qlambda.asDiagonal();
	X = (X * (t * r)).array().exp();
	return U * X * U_;
}

void GTR::updateParams() {
	assert(isValid());
	/* decomposite Q = pi_T * R */
	R.setZero();
	for(Matrix4d::Index j = 0; j < Q.cols(); ++j)
		for(Matrix4d::Index i = 0; i < Q.cols(); ++i)
			if(i != j)
				R(i, j) = Q(i, j) / pi(j);
}

void GTR::updateRate() {
	Q = pi.transpose() * R;
	fixRate();
}

void GTR::updateEigenParams() {
	// Decoposite Q with EigenSover
	EigenSolver<Matrix4d> es(Q);
	if(es.info() != Success) {
		cerr << "Rate Matrix Q cannot be solved for Eigen values" << endl;
		abort();
	}
	Qlambda = es.eigenvalues();
	U = es.eigenvectors();
	U_ = U.inverse();
}

} /* namespace EGriceLab */
