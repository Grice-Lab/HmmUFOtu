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
	double beta; // sequence diversity as 1 / (1 - piA^2 - piC^2 - piG^2 - piT^2)
};

inline Matrix4d F81::Pr(double v) const {
	Matrix4d P;
	double e = ::exp(-beta * v);
	for(Matrix4d::Index i = 0; i < P.rows(); ++i)
		for(Matrix4d::Index j = 0; j < P.cols(); ++j)
			P(i, j) = i == j ? e + pi(j) * (1 - e) : pi(j) * (1 - e);

	return P;
}

} /* namespace EGriceLab */

#endif /* SRC_F81_H_ */
