/*
 * DirichletDensity.cpp
 *
 *  Created on: Jun 29, 2016
 *      Author: zhengqi
 */

#include <boost/math/special_functions/digamma.hpp>
#include <iostream>
#include <cmath>
#include "DirichletDensity.h"

namespace EGriceLab {
namespace Math {

using namespace std;
using namespace Eigen;
using boost::math::digamma;

/* static variable definition */
const double DirichletDensity::DEFAULT_ALPHA = 1;
const double DirichletDensity::DEFAULT_WEIGHT = 0;
const string DirichletDensity::FILE_HEADER = "Dirichlet Density Model";

VectorXd DirichletDensity::meanPostP(const VectorXd& freq) const {
	return (freq + alpha) / (freq.sum() + alpha.sum());
}

VectorXd DirichletDensity::weightGradient(const MatrixXd& data) const {
	int K = getK();
	VectorXd grad(K);
	double alphaSum = alpha.sum();
	MatrixXd::Index M = data.cols();
	RowVectorXd nSum = data.colwise().sum();
	for(int i = 0; i < K; ++i) {
		double S = 0;
		for(MatrixXd::Index t = 0; t < M; ++t) {
			S += digamma(static_cast<double> (data(i, t)) + static_cast<double> (alpha(i)) )
					- digamma(static_cast<double> (nSum(t)) + alphaSum);
		}
		grad(i) = alpha(i) * (M * ( digamma(alphaSum) - digamma(static_cast<double> (alpha(i))) ) + S);
	}
	return grad;
}

double DirichletDensity::trainML(const MatrixXd& data) {
	/* initiate the parameters using moment-matctching */
	momentInit(data);
	double c = cost(data);
	for(int it = 0; maxIter <= 0 || it < maxIter; ++it) { // infinite loop to be terminated within
		/* copy old parameters */
		double cOld = c;
		VectorXd alphaOld(alpha);

		VectorXd wGrad = weightGradient(data);
//		cerr << "wGrad:" << wGrad.transpose() << endl;
		/* update weight and parameters */
		w += eta * wGrad;
		alpha = w.array().exp();
		/* check the new parameters for over-fitting */
		if((alpha.array() == 0).any()) {
			cerr << "Potential over-fitting detected. Please choose another MSA training set" << endl;
			return NAN;
		}
		/* calculate new cost */
		c = cost(data);
		double deltaC = cOld - c;
//		fprintf(stderr, "cOld:%lg c:%lg deltaC:%lg\n", cOld, c, deltaC);

		/* termination check */
		if(alpha.isApprox(alphaOld, absEpsParams + relEpsParams * alphaOld.norm())
				&& deltaC >= 0 && deltaC < absEpsCost + relEpsCost * cOld)
			break;
	}
	setTrainingCost(c);
	return c;
}

double DirichletDensity::lpdf(const VectorXd& freq) const {
	assert(freq.size() == alpha.size());
	int K = getK();
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
	out << FILE_HEADER << endl;
	out << "Training cost: " << getTrainingCost() << endl;
	out << "K: " << getK() << endl;
	out << "alpha:" << endl;
	out << alpha.transpose().format(FULL_FORMAT) << endl;
	return out;
}

void DirichletDensity::momentInit(MatrixXd data) {
	int K = getK();
	int M = data.cols();
	if(M < 2)
		return; /* too few freq to estimate */

	/* Normalize the column sum, so the observed data follows Dirichlet-Multinomial distribution */
	double N = data.colwise().sum().maxCoeff();
	for(int t = 0; t < M; ++t)
		data.col(t) *= N / data.col(t).sum();
	/* calculate the Mean (1st-moment) and Var (2nd-moment) of the observed counts */
	VectorXd dataMean = data.rowwise().mean();
	VectorXd dataVar = (data.colwise() - dataMean).rowwise().squaredNorm() / M;
	/* calculate parameter concentration using E(0) and Var(0) */
	double alphaNorm = 0;
	// try each k
	for(int i = 0; i < K; ++i) {
		alphaNorm = (dataVar(i) - N * dataMean(i) + 1) / (dataMean(i) - 1 / N - dataVar(i));
		if(alphaNorm > 0) // a good estimation
			break;
	}
//	cerr << "alphaNorm:" << alphaNorm << endl;

	if(alphaNorm <= 0) // do not use moment initiate
		return;
	/* calculate parameters */
	alpha = dataMean * alphaNorm / N;
	w = alpha.array().log();
}

istream& DirichletDensity::read(istream& in) {
	string line;
	int K;
	double c;
	std::getline(in, line);
	if(line != FILE_HEADER) {
		in.setstate(ios_base::failbit);
		return in;
	}

	std::getline(in, line);
	sscanf(line.c_str(), "Training cost: %lg", &c); /* read in training cost */

	std::getline(in, line);
	sscanf(line.c_str(), "K: %d", &K); /* Read K */

	/* set fields */
	setTrainingCost(c);
	setK(K);
	alpha.resize(K);
	w.resize(K);

	std::getline(in, line); /* ignore alpha line */
	for(VectorXd::Index i = 0; i < K; ++i)
		in >> alpha(i);
	w = alpha.array().log();
	return in;
}

} /* namespace Math */
} /* namespace EGriceLab */

