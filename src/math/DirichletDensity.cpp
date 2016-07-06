/*
 * DirichletDensity.cpp
 *
 *  Created on: Jun 29, 2016
 *      Author: zhengqi
 */

#include <boost/math/special_functions/digamma.hpp>
#include <iostream>
#include "DirichletDensity.h"

namespace EGriceLab {
namespace Math {

using namespace std;
using boost::math::digamma;

VectorXd DirichletDensity::postP(const VectorXd& freq) const {
	return (freq + alpha) / (freq.sum() + alpha.sum());
}

VectorXd DirichletDensity::expGradient(const MatrixXd& data) const {
	int K = getK();
	VectorXd grad(K);
	double alphaSum = alpha.sum();
	MatrixXd::Index M = data.cols();
	VectorXd nSum = data.colwise().sum();
	for(int k = 0; k < K; ++k) {
		double S = 0;
		for(MatrixXd::Index t = 0; t < M; ++t) {
			S += digamma(static_cast<double> (data(k, t)) + static_cast<double> (alpha(k)) )
					- digamma(static_cast<double> (nSum(t)) + alphaSum);
		}
		grad(k) = alpha(k) * (M * ( digamma(alphaSum) - digamma(static_cast<double> (alpha(k))) ) + S);
	}
	return grad;
}

void DirichletDensity::trainML(const MatrixXd& data, double eta, double epsilonCost, double epsilonParams, int maxIt) {
	double c = cost(data);
	for(int it = 0; maxIt <= 0 || it < maxIt; ++it) { // infinite loop to be terminated within
//		cerr << "alpha: " << alpha.transpose() << endl;
		VectorXd expGrad = expGradient(data);
//		cerr << "wGrad:" << expGrad.transpose() << endl;
		/* update weight and parameters */
		w += eta * expGrad;
		VectorXd alphaOld(alpha);
		alpha = w.array().exp();
		/* check the new parameters for over-fitting */
		if((alpha.array() == 0).any()) {
			cerr << "Potential over-fitting detected. Please choose another MSA training set" << endl;
			abort();
		}
//		cerr << "alpha: " << alpha.transpose() << endl;
		/* calculate new cost */
		double cNew = cost(data);
		double deltaC = ::fabs(cNew - c);
//		fprintf(stderr, "c:%lg cNew:%lg deltaC:%lg\n", c, cNew, deltaC);
		/* termination check */
		if(alpha.isApprox(alphaOld, epsilonParams) && deltaC < epsilonCost)
			break;
		c = cNew;
	}
}

double DirichletDensity::lpdf(const VectorXd& freq) const {
	assert(freq.size() == alpha.size());
	int K = alpha.size();
	/* constant part */
	double freqNorm = freq.sum();
	double alphaNorm = alpha.sum();
	double logC = lgamma(freqNorm + 1) + lgamma(alphaNorm) - lgamma(freqNorm + alphaNorm);
	/* product part */
	double logS = 0;
	for(int i = 0; i < K; ++i) {
		logS += lgamma(static_cast<double> (freq(i)) + static_cast<double> (alpha(i)))
				- lgamma(static_cast<double> (freq(i)) + 1)
				- lgamma(static_cast<double> (alpha(i)));
	}
	return logC + logS;
}

ostream& DirichletDensity::print(ostream& out) const {
	out << "Dirichlet Density Model" << endl;
	out << "K: " << getK() << endl;
	out << "alpha: " << alpha.transpose() << endl;
}

istream& DirichletDensity::read(istream& in) {
	return in;
}

} /* namespace Math */
} /* namespace EGriceLab */

