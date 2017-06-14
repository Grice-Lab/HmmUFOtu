/*
 * LinearAlgebraBasic.h
 *
 *  Created on: Jun 16, 2016
 *      Author: zhengqi
 */

#ifndef SRC_MATH_LINEARALGEBRABASIC_H_
#define SRC_MATH_LINEARALGEBRABASIC_H_

#include <Eigen/Dense>
#include <iostream>
#include <cassert>
#include <cmath>

namespace EGriceLab {

namespace Math {

using Eigen::VectorXd;
using Eigen::MatrixBase;

const double NAT2BIT = 1.0 / ::log(2);

/**
 * Normalize a vector
 */
inline VectorXd normalize(const VectorXd& v, double norm = 0) {
	if(norm == 0)
		norm = v.sum();
	return v / v.sum();
}

/*
 * Calculate vector exp
 */
inline VectorXd exp(const VectorXd& v) {
	return v.array().exp();
}

/*
 * Calculate scaled vector exp by a given factor
 */
inline VectorXd scaleExp(const VectorXd& v, double scale) {
	return (v.array() + scale).exp();
}

/*
 * Calculate scaled vector exp by default method
 */
inline VectorXd scaleExp(const VectorXd& v) {
	return scaleExp(v, -v.maxCoeff());
}

/**
 * Calculate Dirichlet PDF given the parameters
 * @param  alpha Direchlet parameters
 * @param  x observed categorized/multinomial value
 * @return  PDF of overserving the data
 */
inline double dDirichlet(const VectorXd& alpha, const VectorXd& x) {
	std::cerr << "alpha:" << alpha.transpose() << " x:" << x.transpose() << std::endl;
	assert(alpha.size() == x.size());
	VectorXd::Index K = alpha.size();
	double sum = x.sum();
	assert(sum > 0);
	VectorXd xNorm = x / sum;

	/* calculate numerator */
	double logNumer = 0;
	for(VectorXd::Index i = 0; i < K; ++i)
		logNumer += (alpha(i) - 1.0) * ::log(xNorm(i));

	/* calculate denominator */
	double logDenom = 0;
	for(VectorXd::Index i = 0; i < K; ++i)
		logDenom += ::lgamma(static_cast<double> (alpha(i))); /* lgamma is a C99 function */
	logDenom -= ::lgamma(alpha.sum());

	std::cerr << "logNumer:" << logNumer << " logDenom:" << logDenom << std::endl;

	return ::exp(logNumer - logDenom);
}


/**
 * calculate relative entropy between a two distribution
 * zero p probs are ignored
 */
inline double relative_entropy(const VectorXd& p, const VectorXd& q) {
	assert(p.rows() == q.rows());
	VectorXd::Index N = p.rows();
	double ent = 0;
	for(VectorXd::Index i = 0; i < N; ++i)
		if(p(i) > 0)
			ent += p(i) * ::log( static_cast<double> (p(i)) / static_cast<double> (q(i)) );
	return NAT2BIT * ent; /* return entropy in bits */
}

} /* end namespace Math */
} /* end namespace EGriceLab */

#endif /* SRC_MATH_LINEARALGEBRABASIC_H_ */
