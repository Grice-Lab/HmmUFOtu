/*
 * DirichletModel.h
 *
 *  Created on: Jun 16, 2016
 *      Author: zhengqi
 */

#ifndef SRC_MATH_DIRICHLETMODEL_H_
#define SRC_MATH_DIRICHLETMODEL_H_

#include <cassert>
#include <cfloat>
#include <iostream>

namespace EGriceLab {
namespace Math {

using std::istream;
using std::ostream;
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
	 * An abstract method to calculate the posterior probability
	 * given the parameters and an observed frequency
	 */
	virtual VectorXd postP(const VectorXd& freq) const = 0;

	/**
	 * Do a maximum likelihood training of all underlying parameters given a training data,
	 * with M columns each an observed frequency vector, and K rows
	 */
	virtual void trainML(const MatrixXd& data, double eta = DEFAULT_ETA, double epsilon = DEFAULT_EPSILON) = 0;

	/**
	 * Calculate the negative gradient of the Dirichlet concentration parameters
	 * using current parameters and observed data at (unbound) exponential scale
	 * Exponential scale is used because the Dirichlet parameters must be strictly positive.
	 */
	virtual VectorXd expGradient(const MatrixXd& freq) const = 0;

	/**
	 * Calculate the cost of observing one data
	 */
	VectorXd cost(const VectorXd& freq) const;

	/**
	 * Calculate the log postP of observing one data
	 */
	VectorXd logPostP(const VectorXd& freq) const;

private:
	/*
	 * internal methods to support input/output method inheritance
	 * these methods should be implemented in all derived classes, but cannot be called directly
	 */
	virtual ostream& print(ostream& out) const = 0;
	virtual istream& read(istream& in) = 0;

public:
	/* non-member friend functions */
	friend ostream& operator<<(ostream& out, const DirichletModel& dm);
	friend istream& operator>>(istream& in, DirichletModel& dm);

private:
	int K; // number of concentration parameters

public:
	static const double DEFAULT_ETA = 0.05; // default step width relative to the gradient used in ML parameter training
	static const int MIN_K = 2; // minimum number of categories
	static const double DEFAULT_EPSILON = FLT_EPSILON;

};

inline VectorXd DirichletModel::cost(const VectorXd& freq) const {
	return -postP(freq).array().log();
}

inline VectorXd DirichletModel::logPostP(const VectorXd& freq) const {
	return postP(freq).array().log();
}

inline ostream& operator<<(ostream& out, DirichletModel& dm) {
	return dm.print(out);
}

inline istream& operator>>(istream& in, DirichletModel& dm) {
	return dm.read(in);
}

} /* namespace Math */
} /* namespace EGriceLab */

#endif /* SRC_MATH_DIRICHLETMODEL_H_ */
