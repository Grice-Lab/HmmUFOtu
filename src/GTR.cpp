/*
 * GTR.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: zhengqi
 */

#include <cassert>
#include "GTR.h"

namespace EGriceLab {
using namespace std;
using namespace Eigen;

const string name = "GTR";

Matrix4d GTR::Pr(double t, double r) const {
	Matrix4d X = Qlambda.asDiagonal();
	X = (X * (t * r)).array().exp();
	return U * X * U_inv;
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
	U_inv = U.inverse();
}

void GTR::initParam() {
	/* normalization constant */
	double Z = alpha + beta + gamma + delta + epsilon + eta;
	assert(Z > 0);
	alpha /= Z;
	beta /= Z;
	gamma /= Z;
	delta /= Z;
	epsilon /= Z;
	eta /= Z;
	R(A, C) = R(C, A) = alpha;
	R(A, G) = R(G, A) = beta;
	R(A, T) = R(T, A) = gamma;
	R(C, G) = R(G, C) = delta;
	R(C, T) = R(T, C) = epsilon;
	R(G, T) = R(T, G) = eta;
}

} /* namespace EGriceLab */
