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

const string GTR::name = "GTR";


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

void GTR::trainParams(const vector<Matrix4d>& P_vec, const Vector4d& f) {
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

} /* namespace EGriceLab */
