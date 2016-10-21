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
#include <stdexcept>
#include <iostream>
#include <eigen3/Eigen/Dense>

namespace EGriceLab {
namespace Math {

using std::istream;
using std::ostream;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::IOFormat;

class DirichletModel {
public:
	/* default constructor, do nothing */
	DirichletModel(): K(0) { }

	/* construct a Dirichlet model with given categories */
	explicit DirichletModel(int K): K(K) {
		assert(K >= MIN_K);
	}

	/* virtual destructor, do nothing */
	virtual ~DirichletModel() { }

	/* member methods */
	/**
	 * An abstract method to calculate the posterior probabilities of category
	 * given the parameters and an observed frequency
	 */
	virtual VectorXd meanPostP(const VectorXd& freq) const = 0;

	/**
	 * Initiate the Dirichlet parameters using momenth-matching method,
	 * to get a good starting estimate
	 */
	virtual void momentInit(MatrixXd data) = 0;

	/**
	 * Do a maximum likelihood training of all underlying parameters given a training data,
	 * with M columns each an observed frequency vector, and K rows
	 * return NAN if anything went wrong
	 */
	virtual double trainML(const MatrixXd& data,
			int maxIt = MAX_ITERATION, double eta = DEFAULT_ETA,
			double epsilonCost = DEFAULT_REL_EPS_COST,
			double epsilonParams = DEFAULT_REL_EPS_PARAMS) = 0;

	/**
	 * Calculate the logPDF of observing a data using this model
	 */
	virtual double lpdf(const VectorXd& freq) const = 0;

	/**
	 * Calculate the PDF of observing a data using this model
	 */
	virtual double pdf(const VectorXd& freq) const;

	/**
	 * Calculate the cost of observing an entire data
	 */
	double cost(const MatrixXd& data) const;

	/*
	 * internal methods to support input/output method inheritance
	 */
	virtual ostream& print(ostream& out) const = 0;
	virtual istream& read(istream& in) = 0;

public:
	/* non-member friend functions */
	friend istream& operator>>(istream& in, DirichletModel& dm);
	friend ostream& operator<<(ostream& out, const DirichletModel& dm);

	/* getters and setters */
	int getK() const {
		return K;
	}

	/**
	 * Set K # of categories
	 */
	virtual void setK(int k) {
		K = k;
	}

private:
	int K; // number of parameters

	/* static members */
public:
	static const double DEFAULT_ETA = 0.001; // default step width relative to the gradient used in ML parameter training
	static const int MIN_K = 2; // minimum number of categories
//	static const double DEFAULT_EPSILON = FLT_EPSILON;
	static const double DEFAULT_REL_EPS_COST = 1e-6;
	static const double DEFAULT_REL_EPS_PARAMS = 1e-6;
	static const int MAX_ITERATION = 0;
	static const IOFormat FULL_FORMAT; /* default output format for eigen objects */
};

inline ostream& operator<<(ostream& out, const DirichletModel& dm) {
	return dm.print(out);
}

inline istream& operator>>(istream& in, DirichletModel& dm) {
	return dm.read(in);
}

inline double DirichletModel::pdf(const VectorXd& data) const {
	return ::exp(lpdf(data));
}

inline double DirichletModel::cost(const MatrixXd& data) const {
	double c = 0;
	MatrixXd::Index M = data.cols();
	for(MatrixXd::Index t = 0; t < M; ++t) {
		c -= lpdf(data.col(t));
	}
	return c;
}

} /* namespace Math */
} /* namespace EGriceLab */

#endif /* SRC_MATH_DIRICHLETMODEL_H_ */
