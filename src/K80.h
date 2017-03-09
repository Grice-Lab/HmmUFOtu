/*
 * K80.h
 *  K80 DNA Substitution Model
 *  Created on: Mar 8, 2017
 *      Author: zhengqi
 */

#ifndef SRC_K80_H_
#define SRC_K80_H_

#include <cmath>
#include "DNASubModel.h"

namespace EGriceLab {

class K80: public DNASubModel {
public:
	/* Constructors */

	/** default constructor */
	K80() : kappa(1) {
		setBeta();
	}

	/* destructor, do nothing */
	virtual ~K80() { }

	/* member methods */
	virtual string modelType() const {
		return name;
	}

	virtual Vector4d getPi() const {
		return pi;
	}

	/**
	 * get the Prob matrix given branch length and optionally rate factor
	 * @override  the base class pure virtual function
	 */
	virtual Matrix4d Pr(double v) const;

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
	virtual K80* clone() const {
		return new K80(*this);
	}

private:
	/** set beta by kappa and pi */
	void setBeta() {
		beta = 1 / (2 * kappa);
	}

	static const string name;
	static const Vector4d pi;

	double kappa; // Ti/Tv ratio
	double beta;  // rate diversity
};

inline Matrix4d K80::Pr(double v) const {
	Matrix4d P;
	double e = ::exp(-4 * beta * v);
	double eV = ::exp(-2 * (1 + kappa) * beta * v);
	P.diagonal().setConstant((1.0 + e + 2 * eV) / 4);
	P(A,G) = P(G,A) = P(C,T) = P(T,C) = (1.0 + e - 2 * eV) / 4;
	P(A,C) = P(A,T) = P(C,A) = P(C,G) = P(G,C) = P(G,T) = P(T,A) = P(T,G) = (1.0 - e) / 4;

	return P;
}

} /* namespace EGriceLab */

#endif /* SRC_K80_H_ */
