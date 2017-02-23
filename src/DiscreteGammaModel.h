/*
 * DiscreteGammaModel.h
 * A Discrete-Gamma distribution model to capture the rate-heterogeinity among different sites
 *  Created on: Feb 17, 2017
 *      Author: zhengqi
 */

#ifndef SRC_DISCRETEGAMMAMODEL_H_
#define SRC_DISCRETEGAMMAMODEL_H_

#include <string>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <math.h> /* C99 header */
#include "HmmUFOtuConst.h"

using Eigen::VectorXd;
using Eigen::VectorXi;

namespace EGriceLab {

class DiscreteGammaModel {
public:
	/* constructors */

	/**
	 * Default constructor
	 */
	DiscreteGammaModel() : alpha(nan), K(0) { }

	/**
	 * Construct a model with given K and alpha
	 */
	DiscreteGammaModel(int K, double alpha) : K(K), alpha(alpha) {
		b.resize(K + 1);
		r.resize(K);
		setBreaks();
		setRates();
	}

	/* member methods */
	int getK() const {
		return K;
	}

	double getShape() const {
		return alpha;
	}

	void setShape(double alpha) {
		this->alpha = alpha;
		setBreaks();
		setRates();
	}

	double rate(int k) const {
		return r(k);
	}

	/**
	 * load model from an input stream
	 */
	istream& load(istream& in);

	/**
	 * save this model to an output stream
	 */
	ostream& save(ostream& out) const;

	/* static methods */
	/**
	 * use moment-matching traning of observed number of mutations to estimate the shape parameter alpha
	 * given gamma-distrubuted rate, the observed number of changes is negative-binomial distribution
	 */
	static double estimateShape(const VectorXi& X);

	/* private member methods */
private:
	/**
	 * Set the break-points according to current alpha
	 */
	void setBreaks();

	/**
	 * Set the average rates of each category according to alpha and breaks
	 */
	void setRates();

	/* member fields */
private:
	double alpha; // shape parameter (and the scare) of the underlying gamma distribution
	int K; // number of discrete categories
	VectorXd b; // break-points to devide Gamma distribution to equal prob-K categories
	VectorXd r; // average rate of each category
};

} /* namespace EGriceLab */

#endif /* SRC_DISCRETEGAMMAMODEL_H_ */
