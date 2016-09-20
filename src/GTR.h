/*
 * GTR.h
 *  A Generalized Time-Reversible DNA Substitution Model
 *  Created on: Apr 29, 2016
 *      Author: zhengqi
 */

#ifndef GTR_H_
#define GTR_H_
#include <string>
#include <eigen3/Eigen/Dense>
#include "DNASubModel.h"

namespace EGriceLab {

using std::string;
using Eigen::Vector4d;
using Eigen::VectorXd;
using Eigen::Matrix4d;

class GTR : public DNASubModel {
public:
	/* Constructors */

	/* destructor, do nothing */
	virtual ~GTR() { }

	/* member methods */
	virtual const string& modelType() const {
		return "GTR";
	}

	/**
	 * get the Prob matrix given branch length and optionally rate factor
	 * @override  the base class pure virtual function
	 */
	Matrix4d Pr(double t, double r = 1.0) const;

	/*
	 * update other model parameters using current rate and pi data
	 * @override  base class pure virtual method
	 */
	virtual void updateParams();

	/*
	 * update rate matrix using current pi and other parameters
	 * @override  base class pure virtual method
	 */
	virtual void updateRate();

private:
	static const string name;

	/* rate parameters, alpha + beta + gamma + delta + epsilon + eta = 1 */
//	double mu; /* substitution rate per site per unit time */

	Matrix4d R; /* Rate parameters, Q = R * pi */

	Vector4d Qlambda; /* stored eigenvalues of Q */
	Matrix4d U; /* matrix with columns as eigen vectors of Q */
	Matrix4d U_; /* U-1 inverse of U */

	/* private methods */
	/** update lambda, U and U_inv based from Q */
	void updateEigenParams();
};

} /* namespace EGriceLab */

#endif /* SRC_GTR_H_ */
