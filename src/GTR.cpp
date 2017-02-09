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
#include "ProgLog.h"

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
			in >> value; // read in model type
			if(value != modelType()) {
				errorLog << "Unmatched Model Type!" << endl;
				errorLog << "Trying to read in a " << value << " model into a " << modelType() << " object" << endl;
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
		else if(tag == "Q:") { // Q section for human read only
			for(Vector4d::Index i = 0; i <= Q.rows(); ++i)
				std::getline(in, line); /* ignore the entire line */
		}
		else {
			errorLog << "Un-recognized line found in GTR Model file: tag: " << tag << endl << line << endl;
			in.setstate(ios_base::badbit);
			return in;
		}
	}
	setQfromParams();
	return in;
}

ostream& GTR::write(ostream& out) const {
	out << "# DNA Substitution Model" << endl;
	out << "Type: " << modelType() << endl;
	out << "pi: " << pi.transpose().format(FULL_FORMAT) << endl;
	out << "R:" << endl << R.format(FULL_FORMAT) << endl;
	out << "Q:" << endl << Q.format(FULL_FORMAT) << endl;
	return out;
}

void GTR::trainParams(const vector<Matrix4d>& Pv, const Vector4d& f) {
	/* estimate pi using mean f */
	pi = f / f.sum();
//	assert(isValidFreq(pi));
//	cerr << "pi estimated: " << pi.transpose() << endl;

	/* estimate Q from Pv using constrained optimization */
	Q.setZero();
	size_t N = 0;
	for(vector<Matrix4d>::const_iterator it = Pv.begin(); it != Pv.end(); ++it) {
		const Matrix4d& Qv = constrainedQfromP(*it);
		if(isValidRate(Qv)) {
			N++;
			Q += scale(Qv);
		}
	}
	Q /= N;
//	cerr << "Q estimated: " << endl << Q << endl;
//	assert(isValidRate(Q));

	/* Decomposite R from Q, as Qij = Rij * pi(j) */
	for(Matrix4d::Index j = 0; j < R.cols(); ++j)
		R.col(j) = Q.col(j) / pi(j);
	R.diagonal().setZero(); /* set diagonal to zeros */

	/* average R to make it symmetric */
	R += R.transpose().eval();
	R /= 2.0;
	/* reset Q */
	setQfromParams();
}

void GTR::setQfromParams() {
	/* reset Q */
	assert(R.diagonal().sum() == 0);
	for(Matrix4d::Index j = 0; j < R.cols(); ++j)
		Q.col(j) = R.col(j) * pi(j);
	/* setting Q's diagnal elements */
	Q.diagonal() = - Q.rowwise().sum();
	Q = scale(Q); /* re-scale Q */

	/* Eigen-decompsite Q
	 * Since Q is guaranteed to be a diagonizable matrix as Q = U * X * U-1
	 * Eigen values and vectors of Q are guaranteed to be real (not complex)
	 */
	EigenSolver<Matrix4d> es(Q);
	if(es.info() != Eigen::Success) {
		errorLog << "Cannot perform EigenSolver on rate matrix Q:" << endl << Q << endl;
		abort();
	}
	lambda = es.eigenvalues().real();
	U = es.eigenvectors().real();
	U_1 = U.inverse();
}

} /* namespace EGriceLab */
