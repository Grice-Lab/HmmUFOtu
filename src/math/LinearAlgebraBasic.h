/*
 * LinearAlgebraBasic.h
 *
 *  Created on: Jun 16, 2016
 *      Author: zhengqi
 */

#ifndef SRC_MATH_LINEARALGEBRABASIC_H_
#define SRC_MATH_LINEARALGEBRABASIC_H_

#include <eigen3/Eigen/Dense>

namespace EGriceLab {

namespace Math {

using Eigen::VectorXd;

/**
 * Normalize a vector
 */
VectorXd normalize(const VectorXd& v) {
	return v / v.sum();
}

/*
 * Calculate vector exp
 */
VectorXd exp(const VectorXd& v) {
	return v.array().exp();
}

/*
 * Calculate scaled vector exp by a given factor
 */
VectorXd scaleExp(const VectorXd& v, double scale) {
	return (v.array() + scale).exp();
}

/*
 * Calculate scaled vector exp by default method
 */
VectorXd scaleExp(const VectorXd& v) {
	return scaleExp(v, -v.maxCoeff());
}

} /* end namespace Math */
} /* end namespace EGriceLab */

#endif /* SRC_MATH_LINEARALGEBRABASIC_H_ */
