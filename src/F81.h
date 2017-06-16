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
 * F81.h
 *  F81 DNA Substitution Model
 *  Created on: Mar 7, 2017
 *      Author: zhengqi
 */

#ifndef SRC_F81_H_
#define SRC_F81_H_

#include <cmath>
#include "DNASubModel.h"

namespace EGriceLab {

class F81: public DNASubModel {
public:
	/* Constructors */

	/** default constructor */
	F81() : pi(Vector4d::Constant(1.0/4))
	{
		setBeta();
	}

	/* destructor, do nothing */
	virtual ~F81() { }

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
	 * the actual formula is described in McGuire 1999
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
	virtual F81* clone() const {
		return new F81(*this);
	}

private:
	/** set beta by kappa and pi */
	void setBeta() {
		beta = 1 / (1 - pi.squaredNorm());
	}

	static const string name;

	Vector4d pi; /* base frequency */
	double beta; // sequence diversity as 1 / (1 - A^2 - C^2 - G^2 - T^2)
};

inline Matrix4d F81::Pr(double v) const {
	Matrix4d P;
	double e = ::exp(-beta * v);
	for(Matrix4d::Index i = 0; i < P.rows(); ++i)
		for(Matrix4d::Index j = 0; j < P.cols(); ++j)
			P(i, j) = i == j ? e + pi(j) * (1 - e) : pi(j) * (1 - e);

	return P;
}

inline double F81::subDist(const Matrix4d& D, double N) const {
	if(N == 0)
		return 0;
	double p = (D.sum() - D.diagonal().sum()) / N;
	double E = 1 - pi.squaredNorm();
	return - E * ::log(1 - p / E);
}

} /* namespace EGriceLab */

#endif /* SRC_F81_H_ */
