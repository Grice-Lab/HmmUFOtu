/*
 * DirichletMixture.h
 *  Dirichlet Mixture model
 *  Created on: Jul 6, 2016
 *      Author: zhengqi
 */

#ifndef SRC_MATH_DIRICHLETMIXTURE_H_
#define SRC_MATH_DIRICHLETMIXTURE_H_

#include <string>
#include "DirichletModel.h"

namespace EGriceLab {
namespace Math {

using std::string;

class DirichletMixture: public DirichletModel {
public:
	/* constructors */
	/* default constructor */
	DirichletMixture() : L(0) { }

	/* construct a Dirichlet density with given categories and optionally estimated alpha */
	DirichletMixture(int K, int L):
		DirichletModel(K), L(L), q(L), w(K, L), alpha(K, L) /* initiate w and alpha to correct size */ {
		assert(K >= MIN_K);
		assert(L >= MIN_COMPONENT);
		alpha.setConstant(DEFAULT_ALPHA);
		w.setConstant(DEFAULT_WEIGHT);
		q.setConstant(1.0 / L);
	}

	/* destructor, do nothing */
	virtual ~DirichletMixture() { }

	/* member methods */
	/**
	 * Set K
	 * @param K  # of categories
	 * @override  base class virtual method
	 */
	virtual void setK(int K);

	/**
	 * Set L
	 * @param L:  # of mixtures
	 */
	void setL(int L);

	/**
	 * Set dimentions
	 * @param K:  # of categories
	 * @param L:  # of mixtures
	 */
	void setDims(int K, int L);

	/**
	 * Calculate the mean posterior probability given this model an observed frequency
	 * implement the base case abstract method
	 */
	virtual VectorXd meanPostP(const VectorXd& freq) const;

	/**
	 * Calculate the posterior probability of each component given the model parameters and an observed frequency
	 */
	VectorXd compPostP(const VectorXd& freq) const;

	/**
	 * Calculate the negative gradient of the weights (exp(parameters))
	 * using current parameters and observed data
	 * @return  the weight gradient matrix (K x L)
	 */
	MatrixXd weightGradient(const MatrixXd& data) const;

	/**
	 * Initiate the Dirichlet parameters using momenth-matching method
	 * Implement the base class method
	 */
	virtual void momentInit(MatrixXd data);

	/**
	 * Do a maximum likelihood training of all parameters given a training data,
	 * with M columns each an observed frequency vector of length and K (K * M matrix)
	 * @return  cost at trained parameters, or NAN if over-fitting or numeric problem occured
	 */
	virtual double trainML(const MatrixXd& data,
			int maxIt = MAX_ITERATION, double eta = DEFAULT_ETA,
			double epsilonCost = DEFAULT_REL_EPS_COST,
			double epsilonParams = DEFAULT_REL_EPS_PARAMS);

	/**
	 * Calculate the PDF of observing a data using this DM model
	 * the lpdf must base on this method because the existence of the mixture coefficients
	 */
	double pdf(const VectorXd& freq) const;

	/**
	 * Calculate the log PDF of observing a data using this model
	 * implemet the base abstract method
	 */
	virtual double lpdf(const VectorXd& freq) const {
		return ::log(pdf(freq));
	}

	/* static methods */
	/**
	 * log Beta function defined on a vector value
	 * Beta(X) = Pi(Gamma(xi)) / Gamma(|xi|)
	 */
	static double lbeta(const VectorXd& x);

	/* implement base class private method */
	virtual ostream& print(ostream& out) const;
	virtual istream& read(istream& in);

	/* member fields */
private:
	int L; // number of Dirichlet Mixture components
	MatrixXd alpha; // Dirichlet Mixture parameters (K x L), with each column a Dirichlet parameter vector of length K
	MatrixXd w; // Dirichlet Mixture weights (K x L), alpha = exp(w)
	VectorXd q; // mixture coefficient with length L

public:
	static const int MIN_COMPONENT = 2; /* minimum number of components */
	static const double DEFAULT_ALPHA = 1;
	static const double DEFAULT_WEIGHT = 0;
	static const string FILE_HEADER;
};

inline void DirichletMixture::setK(int K) {
	if(K < MIN_K)
		throw std::invalid_argument("DirichletDensity K must be at least " + MIN_K);

	DirichletModel::setK(K); // invoke base class method
	alpha.resize(K, L);
	w.resize(K, L);
	q.resize(L);
	alpha.setConstant(DEFAULT_ALPHA);
	w.setConstant(DEFAULT_WEIGHT);
	q.setConstant(1.0 / L);
}

inline void DirichletMixture::setL(int L) {
	if(L < MIN_COMPONENT)
		throw std::invalid_argument("DirichletMixture L must be at least " + MIN_COMPONENT);

	alpha.resize(getK(), L);
	w.resize(getK(), L);
	q.resize(L);
	alpha.setConstant(DEFAULT_ALPHA);
	w.setConstant(DEFAULT_WEIGHT);
	q.setConstant(1.0 / L);
}

inline void DirichletMixture::setDims(int K, int L) {
	if(K < MIN_K)
		throw std::invalid_argument("DirichletDensity K must be at least " + MIN_K);
	if(L < MIN_COMPONENT)
		throw std::invalid_argument("DirichletMixture L must be at least " + MIN_COMPONENT);

	DirichletModel::setK(K); // invoke base class method
	this->L = L;
	alpha.resize(K, L);
	w.resize(K, L);
	q.resize(L);
	alpha.setConstant(DEFAULT_ALPHA);
	w.setConstant(DEFAULT_WEIGHT);
	q.setConstant(1.0 / L);
}

} /* namespace Math */
} /* namespace EGriceLab */

#endif /* SRC_MATH_DIRICHLETMIXTURE_H_ */
