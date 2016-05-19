/*
 * GTR.h
 *
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

	/** get the Prob matrix given branch length and optionally rate factor
	 */
	Matrix4d Pr(double t, double r = 1.0) const;

	/** prepare rate parameters R */
	virtual void updateParam(const Matrix4d& freq);

private:
	static const string name;

	/* rate parameters, alpha + beta + gamma + delta + epsilon + eta = 1 */
	double mu; /* substitution rate per site per unit time */
	double alpha; /* T_ac = T_ca */
	double beta; /* T_ag = T_ga */
	double gamma; /* T_at = T_ta */
	double delta; /* T_cg = T_gc */
	double epsilon; /* T_ct = T_tc */
	double eta; /* T_gt = T_tg */

	Vector4d Qlambda; /* eigenvalues of Q */
	Matrix4d U; /* matrix with columns as eigen vectors of Q */
	Matrix4d U_inv; /* U-1 inverse of U */

	/* private methods */
	/** update lambda, U and U_inv based from Q */
	void updateEigenParams();
};

} /* namespace EGriceLab */

#endif /* SRC_GTR_H_ */
