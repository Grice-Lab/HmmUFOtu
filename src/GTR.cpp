/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
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
namespace HmmUFOtu {

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
			for(Vector4d::Index i = 0; i < 4; ++i)
				for(Vector4d::Index j = 0; j < 4; ++j)
					in >> R(i, j);
		}
		else if(tag == "Q:") { // Q section for human read only
			for(Vector4d::Index i = 0; i <= Q.rows(); ++i)
				std::getline(in, line); /* ignore the entire line */
			break;
		}
		else {
			errorLog << "Un-recognized line found in GTR Model input: tag: " << tag << endl << line << endl;
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

double GTR::subDist(const Matrix4d& D, double N) const {
	if(N == 0)
		return 0;
	/* get F from D */
	Matrix4d F = D / N;
	Matrix4d Fnorm = (F + F.transpose()) / 2;
	Matrix4d P = pi.asDiagonal() * Fnorm; /* get P using symmetric F */
	/* normalize P by rows */
	P.array().colwise() /= P.rowwise().sum().array();

	/* do matrix log by diagonalizable matrix decomposition */
	EigenSolver<Matrix4d> es(P);
	if(es.info() != Eigen::Success) {
		cerr << "Cannot perform EigenSolver on observed frequency data P:" << endl << P << endl;
		abort();
	}
	Vector4d PSI = es.eigenvalues().real();
	Matrix4d U = es.eigenvectors().real();
	Matrix4d U_1 = U.inverse();

	return - (U * PSI.array().log().matrix().asDiagonal() * U_1).trace();
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
