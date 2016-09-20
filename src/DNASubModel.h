/*
 * DNASubModel.h
 * An abstract class providing interface and basic methods for a DNA Substitution Model
 * static utility methods are also provided
 *  Created on: Apr 1, 2016
 *      Author: zhengqi
 */

#ifndef DNASUBMODEL_H_
#define DNASUBMODEL_H_
#include <string>
#include <iostream>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include "DegenAlphabet.h"
#include "MSA.h"
#include "StringUtils.h"
#include "PhyloTree.h"

namespace EGriceLab {

using std::string;
using std::istream;
using std::ostream;
using Eigen::Vector4d;
using Eigen::VectorXd;
using Eigen::Matrix4d;

class DNASubModel {
public:
	/* Constructors */
	/* virtual destructor, do nothing */
	virtual ~DNASubModel() { }

	/* member methods */
	/* getters and setters */

	/**
	 * set rate parameter matrix from a given matrix
	 */
	void setRate(const Matrix4d& data) {
		assert(isValidRate(data));
		Q = data;
	}

	/* set base frequency vector using observed frequency */
	void setFreq(const Vector4d& freq) {
		pi = freq;
	}

	/** get model type */
	virtual const string& modelType() const = 0;

	/** Scale the rate matrix Q so that a branch length of 1 yields one expected change */
	void scale(double mu = 1.0);

	void fixRate();

	/* update other model parameters using current rate and pi data */
	virtual void updateParams() = 0;

	/* update rate matrix from current parameters */
	virtual void updateRate() = 0;

	/** pure virtual method to get the Prob matrix given branch length and optionally rate factor
	 * @param t  branch length in the unit time
	 * @param r  rate factor in variable rate ML estimates
	 * @return  Probability rate matrix between for Bases
	 */
	virtual Matrix4d Pr(double t, double r) const = 0;

	bool isValid() const {
		return isValidRate(Q) && isValidFreq(pi);
	}

	/* IO methods */
protected:
	/**
	 * read in content from input
	 */
	virtual istream& read(istream& in) = 0;

	/**
	 * write this model to given output
	 */
	virtual ostream& write(ostream& out) const = 0;

public:
	/* static methods */
	/** calculate the p-distance between two aligned DigitalSeq */
	static double pDist(const DigitalSeq& seq1, const DigitalSeq& seq2);

	/** calculate the observed transition frequencies using Goldman (two-sequence) method */
	static Matrix4d calcTransFreq2Seq(const DigitalSeq& seq1, const DigitalSeq& seq2);

	/** Update the parameters using Gojobori (three-sequence) method */
	static Matrix4d calcTransFreq3Seq(const DigitalSeq& outer,
			const DigitalSeq& seq1, const DigitalSeq& seq2);

	/** Scale a rate matrix Q so that a branch length of 1 yields mu expected change in a unit time */
	static Matrix4d scale(Matrix4d Q, Vector4d pi = Vector4d::Ones(), double mu = 1.0);

	/**
	 * Obtain substitution Rate matrix Q from observed frequency matrix using matrix-log method
	 * This method might generate non-valid Q that has negative off-diagnal elements
	 */
	static Matrix4d logQfromP(Matrix4d P, bool reversible = true);

	/**
	 * Obtain substitution Rate matrix Q from observed frequency matrix using matrix-log method
	 * This method might generate non-valid Q that has negative off-diagnal elements
	 */
	static Matrix4d constrainedQfromP(Matrix4d P, bool reversible = true);

	/**
	 * Test whether a 4x4 matrix is a valid rate matrix
	 * A rate matrix requires non-negative off-diagonal elements
	 */
	static bool isValidRate(const Matrix4d& Q);

	/**
	 * Test whether a vector is a valid frequency vector
	 * A valid freq vector must be non-negative and sum to 1
	 */
	static bool isValidFreq(const Vector4d& pi);

	/* friend functions */
	friend istream& operator>>(istream& in, DNASubModel& model);

	friend ostream& operator<<(ostream& out, const DNASubModel& model);

private:
	Vector4d pi; /* base frequency of ATCG */
	Matrix4d Q; /* rate matrix */

public:
	static const DegenAlphabet* abc = SeqCommons::nuclAbc;

};

inline void DNASubModel::setFreq(const Vector4d& freq) {
	pi = freq / freq.sum();
}

inline void DNASubModel::scale(double mu) {
	Q = scale(Q, pi, mu);
}

inline void DNASubModel::fixRate() {
	for(Matrix4d::Index i = 0; i < Q.rows(); ++i) {
		double sum = 0;
		for(Matrix4d::Index j = 0; j < Q.cols(); ++j)
			if(i != j)
				sum += Q(i, j);
		Q(i, i) = - sum;
	}
}


inline bool DNASubModel::isValidRate(const Matrix4d& Q) {
	for(Matrix4d::Index i = 0; i < Q.rows(); ++i)
		for(Matrix4d::Index j = 0; j < Q.cols(); ++j)
			if(i != j && Q(i, j) < 0)
				return false;
	return true;
}

inline bool DNASubModel::isValidFreq(const Vector4d& pi) {
	return (pi.array() >= 0).all() && pi.sum() == 1.0;
}

inline istream& operator>>(istream& in, DNASubModel& model) {
	return model.read(in);
}

inline ostream& operator<<(ostream& out, const DNASubModel& model) {
	return model.write(out);
}

} /* namespace EGriceLab */

#endif /* SRC_DNASUBMODEL_H_ */
