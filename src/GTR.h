/*
 * GTR.h
 *  A Generalized Time-Reversible DNA Substitution Model
 *  Created on: Apr 29, 2016
 *      Author: zhengqi
 */

#ifndef GTR_H_
#define GTR_H_
#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <cassert>
#include "DNASubModel.h"

namespace EGriceLab {

using std::string;
using Eigen::Vector4d;
using Eigen::VectorXd;
using Eigen::Matrix4d;
using Eigen::Vector4cd;
using Eigen::Matrix4cd;

class GTR : public DNASubModel {
public:
	/* Constructors */

	/* destructor, do nothing */
	virtual ~GTR() { }

	/* member methods */
	virtual string modelType() const {
		return "GTR";
	}

	virtual Vector4d getPi() const {
		return pi;
	}

	/**
	 * get the Prob matrix given branch length and optionally rate factor
	 * @override  the base class pure virtual function
	 */
	Matrix4d Pr(double t) const;

	/**
	 * read in content from input stream
	 * will set badbit if anything went wrong
	 * @override  base class method
	 */
	virtual istream& read(istream& in);

	/**
	 * write this model to given output stream
	 * @override  base class method
	 */
	virtual ostream& write(ostream& out) const;

	/**
	 * train model parameters using given sets of observed base transition and frequency counts
	 * @override  base class method
	 */
	virtual void trainParams(const vector<Matrix4d>& Pv, const Vector4d& f);

	/**
	 * copy this object and return the new object's address
	 * @override  base class method
	 */
	virtual GTR* clone() const {
		return new GTR(*this);
	}

private:
	static const string name;

	/* rate parameters, alpha + beta + gamma + delta + epsilon + eta = 1 */
//	double mu; /* substitution rate per site per unit time */

	Vector4d pi; /* base frequency */
	Matrix4d Q; /* Rate matrix */
	Matrix4d R; /* Symmetric rate parameters, Q = pi_T * R for i != j; R(i,i) = 0 */

	Vector4d lambda; /* stored eigenvalues of Q for fast computation */
	Matrix4d U; /* stored eigen-matrix with columns as eigen vectors of Q */
	Matrix4d U_1; /* U-1 inverse of U */

	void setQfromParams();
};

inline Matrix4d GTR::Pr(double t) const {
	assert(t >= 0);
	if(t == 0)
		return Matrix4d::Identity(); /* identity matrix */
	return U * (lambda * t).array().exp().matrix().asDiagonal() * U_1;
}

} /* namespace EGriceLab */

#endif /* GTR_H_ */
