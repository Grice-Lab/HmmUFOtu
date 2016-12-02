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
#include <stdexcept>
#include <eigen3/Eigen/Dense>
#include "DegenAlphabet.h"
#include "MSA.h"
#include "StringUtils.h"
#include "PhyloTree.h"
#include "HmmUFOtuConst.h"

namespace EGriceLab {

using std::string;
using std::istream;
using std::ostream;
using std::vector;
using Eigen::Vector4d;
using Eigen::VectorXd;
using Eigen::Matrix4d;
using Eigen::IOFormat;

class DNASubModel {
public:
	/* Constructors */
	/* virtual destructor, do nothing */
	virtual ~DNASubModel() { }

	/* member methods */
	/** get model type */
	virtual string modelType() const = 0;

	/** Get the P transition matrix given a branch length (in unit time ) and an optional rate factor
	 * @param t  branch length in the unit time
	 * @param r  rate factor in variable rate ML estimates
	 * @return  Probability rate matrix between for Bases
	 */
	virtual Matrix4d Pr(double t, double r = 1.0) const = 0;

	/**
	 * train model parameters using a given Phylogenetic tree and method
	 */
	void trainParams(const PhyloTree& tree, string method = "Gojobori");

	/**
	 * Evaluate the likelihood (cost) a Phylogenetic tree using this model
	 * All tree sites are assumed to be invariated
	 * @param tree  tree to evaluate, interval values will be set
	 */
	void evaluate(PhyloTree& tree) const;

	/**
	 * Evaluate a likelihood (cost) of a Phylogenetic tree at given site using this model
	 * @param tree  tree to evaluate, interval values will be set
	 * @param j  site to evalaute (recursively)
	 */
	virtual void evaluate(PhyloTree& tree, int j) const = 0;

	/** get the total cost of a Phylogenetic tree using this model
	 * Assumes all nodes have been evaluated
	 */
	double cost(const PhyloTree& tree) const;

	/** get the total cost at given position of a Phylogenetic tree using this model
	 * Assumes all nodes have been evaluated
	 */
	virtual double cost(const PhyloTree& tree, int j) const = 0;

protected:
	/* public non-virtual methods that call private virtual methods */
	/**
	 * train model parameters using a given Phylogenetic tree and the "Goldman" method
	 */
	void trainParamsGoldman(const PhyloTree& tree);

	/**
	 * train model parameters using a given Phylogenetic tree and method
	 */
	void trainParamsGojobori(const PhyloTree& tree);

private:
	/**
	 * train model parameters using given sets of observed base transition and overall frequency stored in vector
	 */
	virtual void trainParamsByDataset(const vector<Matrix4d>& P_vec, const Vector4d& f) = 0;

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

	/** calculate the observed transition frequencies using Gojobori (three-sequence) method */
	static Matrix4d calcTransFreq3Seq(const DigitalSeq& outer,
			const DigitalSeq& seq1, const DigitalSeq& seq2);

	/** calculate the observed base frequencies of a given seq */
	static Vector4d calcBaseFreq(const DigitalSeq& seq);

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
	static bool isValidRate(Matrix4d Q);

	/**
	 * Test whether a vector is a valid frequency vector
	 * A valid freq vector must be non-negative and sum to 1
	 */
	static bool isValidFreq(const Vector4d& pi);

	/* friend functions */
	friend istream& operator>>(istream& in, DNASubModel& model);

	friend ostream& operator<<(ostream& out, const DNASubModel& model);

public:
	static const double MAX_PDIST; /* maximum p-dist between training sequences */
	static const IOFormat FULL_FORMAT; /* default output format for eigen objects */
	static const IOFormat STD_FORMAT; /* standard output format for eigen objects */
};

inline void DNASubModel::trainParams(const PhyloTree& tree, string method) {
	if(StringUtils::toLower(method) == "goldman")
		return trainParamsGoldman(tree);
	else if(StringUtils::toLower(method) == "gojobori")
		return trainParamsGojobori(tree);
	else
		throw std::invalid_argument("Unknown DNA model training method '" + method + "'");
}

inline void DNASubModel::evaluate(PhyloTree& tree) const {
	for(int j = 0; j < tree.alnSites(); ++j)
		evaluate(tree, j);
}

inline double DNASubModel::cost(const PhyloTree& tree) const {
	double c = 0;
	for(int j = 0; j < tree.alnSites(); ++j)
		c += cost(tree, j);
	return c;
}

inline bool DNASubModel::isValidRate(Matrix4d Q) {
	/* set the diagonal to zeros of this copy */
	if((Q.array() == 0).all()) /* all zero rate is invalid */
		return false;
	Q.diagonal().setZero();
	return (Q.array() >= 0).all();
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

#endif /* DNASUBMODEL_H_ */
