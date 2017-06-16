/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
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
	 * Get the substitution distance given the observed fraction of differences (p-distance) using this model
	 * @override  the base class function
	 */
	virtual double subDist(const Matrix4d& D, double N) const;

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

inline double K80::subDist(const Matrix4d& D, double N) const {
	if(N == 0)
		return 0;
	double p = (D(A,G) + D(G,A) + D(C,T) + D(T,C)) / N; /* observed Ti diff */
	double q = (D(A,C) + D(A,T) + D(C,A) + D(C,G) + D(G,C) + D(G,T) + D(T,A) + D(T,G)) / N; /* observed Tv diff */
	return - 1.0 / 2 * ::log(1 - 2 * p - q) - 1.0 / 4 * ::log(1 - 2 * q);
}

} /* namespace EGriceLab */

#endif /* SRC_K80_H_ */
