/*
 * DirichletDensity.h
 *  Single Dirichlet Density model which is a special case for Dirichlet Mixture Model
 *  Created on: Jun 29, 2016
 *      Author: zhengqi
 */

#ifndef SRC_DIRICHLETDENSITY_H_
#define SRC_DIRICHLETDENSITY_H_

#include "DirichletModel.h"
#include "LinearAlgebraBasic.h"

namespace EGriceLab {
namespace Math {

class DirichletDensity: public DirichletModel {
public:
	/* constructors */
	/* default constructor disabled */
/*	DirichletDensity();*/

	/* construct a Dirichlet density with given categories and optionally estimated alpha */
	explicit DirichletDensity(int K, double alphaEst = DEFAULT_ALPHA):
		DirichletModel(K), w(K), alpha(K) /* initiate w and alpha to correct size */ {
		assert(K >= MIN_K);
		assert(alphaEst > 0);
		alpha.setConstant(alphaEst);
		w = alpha.array().log();
	}

	/* destructor, do nothing */
	virtual ~DirichletDensity() { }

	/* member methods */
	/**
	 * Calculate the posterior probability given this model an observed frequency
	 * implement the base case abstract method
	 */
	virtual VectorXd postP(const VectorXd& freq) const;

	/**
	 * Calculate the negative gradient of the weights (exp(parameters))
	 * using current parameters and observed data
	 */
	virtual VectorXd expGradient(const MatrixXd& data) const;

	/**
	 * Do a maximum likelihood training of all parameters given a training data,
	 * with M columns each an observed frequency vector of length and K (K * M matrix)
	 */
	virtual void trainML(const MatrixXd& data,
			double eta = DEFAULT_ETA, double epsilonCost = DEFAULT_EPSILON_COST,
			double epsilonParams = DEFAULT_EPSILON_PARAMS, int maxIt = MAX_ITERATION);

	/**
	 * Calculate the log PDF of observing a data using this model
	 * implemet the base abstract method
	 */
	virtual double lpdf(const VectorXd& freq) const;

private:
	/* implement base class private method */
	virtual ostream& print(ostream& out) const;
	virtual istream& read(istream& in);

	/* member fields */
private:
	VectorXd alpha; // Dirichlet density parameters
	VectorXd w; // weight parameters, alpha = exp(w)

public:
	static const double DEFAULT_ALPHA = 1;
};

} /* namespace Math */
} /* namespace EGriceLab */

#endif /* SRC_DIRICHLETDENSITY_H_ */
