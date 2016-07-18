/*
 * DNASubModel.h
 *
 *  Created on: Apr 1, 2016
 *      Author: zhengqi
 */

#ifndef DNASUBMODEL_H_
#define DNASUBMODEL_H_
#include <string>
#include <eigen3/Eigen/Dense>
#include "DegenAlphabet.h"
#include "MSA.h"
#include "StringUtils.h"
#include "PhyloTree.h"

namespace EGriceLab {

using std::string;
using Eigen::Vector4d;
using Eigen::VectorXd;
using Eigen::Matrix4d;

class DNASubModel {
public:
	/* nested class and enums */
	enum Base { A, C, G, T };

public:
	/* Constructors */
	/* virtual destructor, do nothing */
	virtual ~DNASubModel() { }

	/* member methods */
	/** get model name */
	string getName() const;

	/** pure virtual method to prepare rate parameters R */
	virtual void updateParam(const Matrix4d& freq) = 0;

	/* init (normalize) base frequency vector pi */
	void initFreq();

	/** update Q from R and pi */
	void updateRate();

	/** Scale the Q so that a branch lenth of 1 yields one expected change */
	void scaleRate();

	/** Estimate the substitution frequency from an multiple sequence alignment using a given algorithm */
	Matrix4d estimateSubRate(const PhyloTree* tree, const string& method = "Gojobori");

protected:
	/** Estimate the substitution parameters using Goldman algorithm */
	Matrix4d estimateSubRateGoldman(const PhyloTree* tree);

	/** Estimate the substitution parameters using Goldman algorithm */
	Matrix4d estimateSubRateGojobori(const PhyloTree* tree);

	/** Update the parameters using two-sequence method */
	virtual Matrix4d updateParams2Seq(const PhyloTree::PhyloTreeNode* seq1, const PhyloTree::PhyloTreeNode* seq2);
	/** Update the parameters using three-sequence method */
	virtual Matrix4d updateParams3Seq(const PhyloTree::PhyloTreeNode* outer,
			const PhyloTree::PhyloTreeNode* seq1, const PhyloTree::PhyloTreeNode* seq2);

public:
	/** pure virtual method to get the Prob matrix given branch length and optionally rate factor
	 * @param v  branch length in the unit time
	 * @param h  rate factor in variable rate ML estimates
	 * @return  Probability rate matrix between for Bases
	 */
	virtual Matrix4d Pr(double t, double r) const = 0;


private:
	string name;

	Vector4d pi; /* base frequency of ATCG */
	Matrix4d R; /* Symmetric rate parameter matrix */
	Matrix4d Q; /* rate matrix, Q = pi X R for i!=j */
};

inline string DNASubModel::getName() const {
	return name;
}

inline void DNASubModel::initFreq() {
	pi /= pi.sum();
}

inline Matrix4d DNASubModel::estimateSubRate(const PhyloTree *tree, const string& method) {
	if(StringUtils::toLower(method) == "goldman")
		return estimateSubRateGoldman(tree);
	else if(StringUtils::toLower(method) == "gojobori")
		return estimateSubRateGojobori(tree);
	else {
		std::cerr << "Unknown parameter estimation method " << method << std::endl;
		abort();
		return Matrix4d::Zero();
	}
}

} /* namespace EGriceLab */

#endif /* SRC_DNASUBMODEL_H_ */
