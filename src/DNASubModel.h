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

namespace EGriceLab {

using std::string;
using Eigen::Vector4d;
using Eigen::VectorXd;
using Eigen::Matrix4d;

class DNASubModel {
public:
	/* Constructors */
	/* virtual destructor, do nothing */
	virtual ~DNASubModel() { }

	/* member methods */
	/** get model name */
	string getName() const;

	/**
	 * pure virtual method to prepare rate parameters R
	 * @param freq  observed transition frequencies between different bases
	 */
	virtual void updateParam(const Matrix4d& freq) = 0;

	/* init (normalize) base frequency vector pi */
	void initFreq();

	/** update Q from R and pi */
	void updateRate();

	/** Scale the Q so that a branch lenth of 1 yields one expected change */
	void scaleRate();

	/** Estimate the substitution frequency from an multiple sequence alignment using a given algorithm */
//	Matrix4d estimateSubRate(const PhyloTree* tree, const string& method = "Gojobori");

	/** Update the parameters using two-sequence method */
	virtual Matrix4d estimateParams2Seq(const DigitalSeq& seq1, const DigitalSeq& seq2);
	/** Update the parameters using three-sequence method */
	virtual Matrix4d estimateParams3Seq(const DigitalSeq& outer,
			const DigitalSeq& seq1, const DigitalSeq& seq2);

public:
	/** pure virtual method to get the Prob matrix given branch length and optionally rate factor
	 * @param v  branch length in the unit time
	 * @param h  rate factor in variable rate ML estimates
	 * @return  Probability rate matrix between for Bases
	 */
	virtual Matrix4d Pr(double t, double r) const = 0;

	/* static DNA alphabet used to encode/decode all DNA Sub models */
	static const DegenAlphabet* abc = SeqCommons::nuclAbc;

private:
	string name;

	Vector4d pi; /* base frequency of ATCG */
	Matrix4d R; /* rate parameters */
	Matrix4d Q; /* rate matrix, Q = pi X R for i!=j */

	/* static methods */
public:
	/**
	 * Obtain substitution Rate matrix Q from observed frequency matrix using matrix-log method
	 * This method might generate non-valid Q that has negative off-diagnal elements
	 */
	static Matrix4d logQfromP(const Matrix4d P);

	/**
	 * Obtain substitution Rate matrix Q from observed frequency matrix using matrix-log method
	 * This method might generate non-valid Q that has negative off-diagnal elements
	 */
	static Matrix4d constrainedQfromP(const Matrix4d P);
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
