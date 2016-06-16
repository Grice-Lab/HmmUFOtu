/*
 * DirichletModel.h
 *
 *  Created on: Jun 16, 2016
 *      Author: zhengqi
 */

#ifndef SRC_MATH_DIRICHLETMODEL_H_
#define SRC_MATH_DIRICHLETMODEL_H_

#include <cassert>

namespace EGriceLab {
namespace Math {

using Eigen::VectorXd;
using Eigen::MatrixXd;

class DirichletModel {
public:
	/* default constructor, set K to minimum requirement */
	DirichletModel(): K(MIN_K) { }

	/* construct a Dirichlet model with given categories */
	explicit DirichletModel(int K): K(K) {
		assert(K >= MIN_K);
	}

	/* virtual destructor, do nothing */
	virtual ~DirichletModel() { }

	/* member methods */
	/**
	 * An abstract method to calculate the posterior probability of an observed multinomial frequency
	 */
	virtual VectorXd postP(const VectorXd& freq) = 0;

	/**
	 * Do a maximum likelihood training of all underlying parameters given a training data,
	 * with m columns each an observed frequency vector, and K rows
	 */
	virtual void trainML(const MatrixXd& data, double eta) = 0;

private:
	int K; // number of concentration parameters

	static const double DEFAULT_ETA = 0.05; // default step width relative to the gradient used in ML parameter training
	static const int MIN_K = 2; // minimum number of categories
};

} /* namespace Math */
} /* namespace EGriceLab */

#endif /* SRC_MATH_DIRICHLETMODEL_H_ */
