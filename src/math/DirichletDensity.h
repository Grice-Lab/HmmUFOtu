/*
 * DirichletDensity.h
 *  Single Dirichlet Density model which is a special case for Dirichlet Mixture Model
 *  Created on: Jun 29, 2016
 *      Author: zhengqi
 */

#ifndef SRC_DIRICHLETDENSITY_H_
#define SRC_DIRICHLETDENSITY_H_

#include <string>
#include <stdexcept>
#include "DirichletModel.h"

namespace EGriceLab {
namespace Math {

using std::string;

class DirichletDensity: public DirichletModel {
public:
	/* constructors */
	/* default constructor, do nothing */
	DirichletDensity() { }

	/* construct a Dirichlet density with given categories and optionally estimated alpha */
	explicit DirichletDensity(int K):
		DirichletModel(K), w(K), alpha(K) /* initiate w and alpha to correct size */ {
		assert(K >= MIN_K);
		alpha.setConstant(DEFAULT_ALPHA);
		w.setConstant(DEFAULT_WEIGHT);
	}

	/* destructor, do nothing */
	virtual ~DirichletDensity() { }

	/* member methods */
	/**
	 * Set K
	 * @param K  # of categories
	 * @override  base class virtual method
	 */
	virtual void setK(int K);

	/**
	 * Calculate the posterior probability given this model an observed frequency
	 * implement the base case abstract method
	 */
	virtual VectorXd meanPostP(const VectorXd& freq) const;

	/**
	 * Calculate the negative gradient of the weights (exp(parameters))
	 * using current parameters and observed data
	 */
	VectorXd weightGradient(const MatrixXd& data) const;

	/**
	 * Initiate the Dirichlet parameters using momenth-matching method
	 * Implement the base class method
	 */
	virtual void momentInit(MatrixXd data);

	/**
	 * Do a maximum likelihood training of all parameters given a training data,
	 * with M columns each an observed frequency vector of length and K (K * M matrix)
	 * implment the base class method
	 * @return  cost at trained parameters, or NAN if anything went wrong
	 */
	virtual double trainML(const MatrixXd& data,
			int maxIt = MAX_ITERATION, double eta = DEFAULT_ETA,
			double epsilonCost = DEFAULT_REL_EPS_COST,
			double epsilonParams = DEFAULT_REL_EPS_PARAMS);

	/**
	 * Calculate the log PDF of observing a data using this model
	 * implemet the base abstract method
	 */
	virtual double lpdf(const VectorXd& freq) const;

	/* implement base class private method */
	virtual ostream& print(ostream& out) const;
	virtual istream& read(istream& in);

	/* member fields */
private:
	VectorXd alpha; // Dirichlet density parameters
	VectorXd w; // weight parameters, alpha = exp(w)

public:
	static const double DEFAULT_ALPHA = 1;
	static const double DEFAULT_WEIGHT = 0;
	static const string FILE_HEADER;
};

inline void DirichletDensity::setK(int K) {
	if(K < MIN_K)
		throw std::invalid_argument("DirichletDensity K must be at least " + MIN_K);
	DirichletModel::setK(K); // invoke base class method
	alpha.resize(K);
	w.resize(K);
	alpha.setConstant(DEFAULT_ALPHA);
	w.setConstant(DEFAULT_WEIGHT);
}

} /* namespace Math */
} /* namespace EGriceLab */

#endif /* SRC_DIRICHLETDENSITY_H_ */
