/*
 * DirichletMixture.cpp
 *
 *  Created on: Jul 6, 2016
 *      Author: zhengqi
 */

#include <cassert>
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
#include <iostream>
#include <algorithm>
#include "DirichletMixture.h"

namespace EGriceLab {
namespace Math {

using namespace std;
using namespace Eigen;
using boost::math::digamma;

/* static variable definition */
const string DirichletMixture::FILE_HEADER = "Dirichlet Mixture Model";

/* private comparator functions */
struct MyFreqComparator {
	MyFreqComparator(const MatrixXd& data) : data(data) { }

	bool operator() (int i, int j) {
		VectorXd f1 = data.col(i);
		VectorXd f2 = data.col(j);
		for(MatrixXd::Index k = 0; k < data.rows(); ++k) {
			if(f1(k) != f2(k))
				return f1(k) < f2(k); // higher priority satisfied
		}
		return false; // all frequencies are equal
	}
	const MatrixXd& data;

};

VectorXd DirichletMixture::meanPostP(const VectorXd& data) const {
	assert(data.size() == alpha.rows());
	int K = alpha.rows();
	/* calculate the beta function part */
	VectorXd logB(L);
	for(int j = 0; j < L; ++j)
		logB(j) = lbeta(alpha.col(j) + data) - lbeta(alpha.col(j));
	logB.array() -= logB.maxCoeff(); // normalize to avoid numeric overflow

	VectorXd X = VectorXd::Zero(K);
	double dataSum = data.sum();
	RowVectorXd alphaSum = alpha.colwise().sum();
	for(int i = 0; i < K; ++i)
		for(int j = 0; j < L; ++j)
			X(i) += q(j) * ::exp(logB(j)) * (alpha(i, j) + data(i)) / (alphaSum(j) + dataSum);
	return X / X.sum();
}

MatrixXd DirichletMixture::weightGradient(const MatrixXd& data) const {
	assert(data.rows() == alpha.rows());
	int K = alpha.rows();
	MatrixXd::Index M = data.cols();

	MatrixXd grad(K, L);
	RowVectorXd alphaSum = alpha.colwise().sum();
	RowVectorXd nSum = data.colwise().sum();

	/* calculate the compPostP on each training data column */
	MatrixXd compP(L, M);
	for(MatrixXd::Index t = 0; t < M; ++t)
		compP.col(t) = compPostP(data.col(t));
	VectorXd compS = compP.rowwise().sum();

	for(int j = 0; j < L; ++j) { // for each component
		for(int i = 0; i < K; ++i) { // for each category
			double S = 0;
			for(MatrixXd::Index t = 0; t < M; ++t) {
				S += compP(j, t) * (digamma(static_cast<double> (data(i, t)) + static_cast<double> (alpha(i, j)))
						- digamma(static_cast<double> (nSum(t)) + static_cast<double> (alphaSum(j)) ));
			}
			grad(i, j) = alpha(i, j) * (compS(j) * (digamma(static_cast<double> (alphaSum(j)))
					- digamma(static_cast<double> (alpha(i, j))) ) + S);
		}
	}
	return grad;
}

double DirichletMixture::trainML(const MatrixXd& data, int maxIt, double eta,
		double epsilonCost, double epsilonParams) {
	assert(data.rows() == alpha.rows());
	/* initiate the parameters using moment-matctching */
	momentInit(data);

	MatrixXd::Index M = data.cols();
	/* EM algorithm to update both the Dirichlet parameters and mixture coefficients */
	double c = cost(data);
	for(int it = 0; maxIt <= 0 || it < maxIt; ++it) { // infinite loop to be terminated within
		/* M step, maximize the parameters using gradient descent */
		MatrixXd wGrad = weightGradient(data);
		/* update weight and parameters */
		w += eta * wGrad;
		MatrixXd alphaOld(alpha);
		alpha = w.array().exp();

		/* check the new parameters for over-fitting */
		if((alpha.array() == 0).any()) {
			cerr << "Potential over-fitting detected. Please choose another MSA training set" << endl;
			return NAN;
		}
		if((q.array() == 0).any()) {
			cerr << "Potential unused (zero-coefficient) mixture component detected. Please use a smaller q for training" << endl;
			return NAN;
		}
		/* calculate new cost */
		double cNew = cost(data);
		double deltaC = (c - cNew) / c;
//		fprintf(stderr, "c:%lg cNew:%lg deltaC:%lg\n", c, cNew, deltaC);

		/* E step, update the mixture coefficients using iteration */
		VectorXd qNew = VectorXd::Zero(L);
		for(MatrixXd::Index t = 0; t < M; ++t)
			qNew += compPostP(data.col(t));
		qNew /= static_cast<double> (M);
		c = cNew;
		q = qNew;

		double alphaNorm = alphaOld.norm();
		/* termination check */
		if(alpha.isApprox(alphaOld, epsilonParams * alphaNorm) && deltaC >= 0 && deltaC < epsilonCost)
			break;
	}
	return c;
}

double DirichletMixture::pdf(const VectorXd& data) const {
	assert(data.size() == alpha.rows());
	int K = alpha.rows();
	double dataNorm = data.sum();

	double p = 0;
	for(int j = 0; j < L; ++j) {
		/* constant part */
		double alphaNorm = alpha.col(j).sum();
		double logC = lgamma(dataNorm + 1) + lgamma(alphaNorm) - lgamma(dataNorm + alphaNorm);
		/* product part */
		double logS = 0;
		for(int i = 0; i < K; ++i) {
			logS += lgamma(static_cast<double> (data(i)) + static_cast<double> (alpha(i, j)))
						- lgamma(static_cast<double> (data(i)) + 1)
						- lgamma(static_cast<double> (alpha(i, j)));
		}
		p += q(j) * ::exp(logC + logS);
	}
	return p;
}

double DirichletMixture::lbeta(const VectorXd& x) {
	double s = 0;
	for(VectorXd::Index i = 0; i != x.size(); ++i)
		s += lgamma(static_cast<double> (x(i)));
	return s - lgamma(x.sum());
}

VectorXd DirichletMixture::compPostP(const VectorXd& data) const {
	assert(data.size() == alpha.rows());
	int K = alpha.rows();
	VectorXd logP(L); // un-normalized component posterior probability
	double dataSum = data.sum();
	RowVectorXd alphaSum = alpha.colwise().sum();
	for(int j = 0; j < L; ++j) {
		double C = lgamma(dataSum + 1) + lgamma(static_cast<double> (alphaSum(j)))
				- lgamma(dataSum + static_cast<double> (alphaSum(j))); // const part
		double S = 0; // product part
		for(int i = 0; i < K; ++i)
			S += lgamma(static_cast<double> (data(i)) + static_cast<double> (alpha(i, j)))
				- lgamma(static_cast<double> (data(i)) + 1)
				- lgamma(static_cast<double> (alpha(i, j)));
		logP(j) = C + S;
	}
	VectorXd p = q.array() * logP.array().exp();
	return p / p.sum();
}

ostream& DirichletMixture::print(ostream& out) const {
	out << FILE_HEADER << endl;
	out << "K: " << getK() << " L:" << L << endl;
	out << "Mixture coefficients:" << endl;
	out << q.transpose().format(FULL_FORMAT) << endl;
	out << "alpha:" << endl;
	out << alpha.format(FULL_FORMAT) << endl;
	return out;
}

void DirichletMixture::momentInit(MatrixXd data) {
	int K = getK();
	int M = data.cols();
	if(M < 2 * L)
		return; /* at least 2 data required for each component */

	/* Sort data columns randomly */
	int* idx = new int[M];
	for(int t = 0; t < M; ++t)
		idx[t] = t;
	std::random_shuffle(idx, idx + M);

	MatrixXd dataSorted(K, M);
	for(int t = 0; t < M; ++t)
		dataSorted.col(t) = data.col(idx[t]);
	delete[] idx;

	/* Normalize the column sum, so the observed data follows Dirichlet-Multinomial distribution */
	double N = dataSorted.colwise().sum().maxCoeff();
	for(int t = 0; t < M; ++t)
		dataSorted.col(t) *= N / dataSorted.col(t).sum();

	/* Divide the data to L equal size categories and do moment-matching */
	for(int j = 0; j < L; ++j) {
		int blockStart = j * M / L;
		MatrixXd block = dataSorted.block(0, blockStart, K, M / L);
		//	cerr << "Calculating moments" << endl;
		/* calculate the Mean (1st-moment) and Var (2nd-moment) of the observed frequencies */
		VectorXd blockMean = block.rowwise().mean();
		VectorXd blockVar = (block.colwise() - blockMean).rowwise().squaredNorm() /block.cols();

		/* calculate parameter concentration try tring E(i) and Var(i) */
		double alphaNorm = 0;
		for(int i = 0; i < K; ++i) {
			alphaNorm = (blockVar(i) - N * blockMean(i) + 1) / (blockMean(i) - 1 / N - blockVar(i));
			if(alphaNorm > 0)
				break;
		}
		if(alphaNorm <= 0) // do not use moment initiate for this component
			continue;

		alpha.col(j) = blockMean * alphaNorm / N;
	}
	w = alpha.array().log();
}

istream& DirichletMixture::read(istream& in) {
	string line;
	int K;
	std::getline(in, line);
	if(line != FILE_HEADER) {
		in.setstate(ios_base::failbit);
		return in;
	}
	std::getline(in, line);
	sscanf(line.c_str(), "K: %d L: %d", &K, &L); /* Read K */
	/* set fields */
	setK(K);
	q.resize(L);
	alpha.resize(K, L);
	w.resize(K, L);

	std::getline(in, line); /* ignore mixture coefficients line */
	for(VectorXd::Index j = 0; j < L; ++j) /* read q */
		in >> q(j);
	std::getline(in, line); /* ignore alpha line */
	for(MatrixXd::Index i = 0; i < K; ++i)
		for(MatrixXd::Index j = 0; j < L; ++j)
			in >> alpha(i, j);
	w = alpha.array().log();
	return in;
}


} /* namespace Math */
} /* namespace EGriceLab */


