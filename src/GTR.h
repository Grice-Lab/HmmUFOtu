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

	/**
	 * get the Prob matrix given branch length and optionally rate factor
	 * @override  the base class pure virtual function
	 */
	Matrix4d Pr(double t, double r = 1.0) const;

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
	 * Evaluate a likelihood (cost) of a Phylogenetic tree at given site using this model
	 * @param tree  tree to evaluate, interval values will be set
	 * @param j  site to evalaute (recursively)
	 * @override  the base class method
	 */
	virtual void evaluate(PhyloTree& tree, int j) const;

	/** get the total cost at given position of a Phylogenetic tree using this model
	 * Assumes all nodes have been evaluated
	 * @override  the base class method
	 */
	virtual double cost(const PhyloTree& tree, int j) const;

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

	/**
	 * train model parameters using given sets of observed base transition and frequency counts
	 * @override  base class method
	 */
	virtual void trainParamsByDataset(const vector<Matrix4d>& P_vec, const Vector4d& f);

	void setQfromParams();
};

inline Matrix4d GTR::Pr(double t, double r) const {
	return U * (lambda * (t * r)).array().exp().matrix().asDiagonal() * U_1;
}

inline double GTR::cost(const PhyloTree& tree, int j) const {
	return ::log(pi.dot(tree.cost.col(j).array().exp().matrix())); /* Eigen guarantee exp(-inf) == 0 */
}

} /* namespace EGriceLab */

#endif /* GTR_H_ */
