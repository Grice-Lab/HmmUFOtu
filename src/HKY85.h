/*
 * HKY85.h
 *
 *  Created on: Mar 7, 2017
 *      Author: zhengqi
 */

#ifndef SRC_HKY85_H_
#define SRC_HKY85_H_

#include <cmath>
#include "DNASubModel.h"

namespace EGriceLab {

class HKY85: public DNASubModel {
public:
	/* Constructors */

	/** default constructor */
	HKY85() : kappa(1), pi(Vector4d::Constant(1.0/4))
	{
		setBeta();
	}

	/* destructor, do nothing */
	virtual ~HKY85() { }

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
	virtual HKY85* clone() const {
		return new HKY85(*this);
	}

private:
	/** set beta by kappa and pi */
	void setBeta() {
		beta = 1 / (2 * (pi(A) + pi(G)) * (pi(C) + pi(T)) + 2 * kappa * (pi(A) * pi(G) + pi(C) * pi(T)));
	}

	static const string name;

	Vector4d pi; /* base frequency */
	double kappa; // Ti/Tv ratio
	double beta; // sequence diversity as 1 / (2(piA + piG)(piC + piT) + 2kata(piA*piG + piC*piT))
};

inline Matrix4d HKY85::Pr(double v) const {
	Matrix4d P;
	double a = pi(A);
	double c = pi(C);
	double g = pi(G);
	double t = pi(T);
	double e = ::exp(-beta * v);
	double eAG = ::exp(-(1 + (a + g) * (kappa - 1)) * beta * v);
	double eCT = ::exp(-(1 + (c + t) * (kappa - 1)) * beta * v);

	P(A, A) = (a * (a + g + (c + t) * e) + g * eAG) / (a + g); /* self */
	P(A, C) = c * (1 - e);                                     /* Tv */
	P(A, G) = (g * (a + g + (c + t) * e) - g * eAG) / (a + g); /* Ti */
	P(A, T) = t * (1 - e);                                     /* Tv */

	P(C, A) = a * (1 - e);                                     /* Tv */
	P(C, C) = (c * (c + t + (a + g) * e) + t * eCT) / (c + t); /* self */
	P(C, G) = g * (1 - e);                                     /* Tv */
	P(C, T) = (t * (c + t + (a + g) * e) - t * eCT) / (c + t); /* Ti */

	P(G, A) = (a * (a + g + (c + t) * e) - a * eAG) / (a + g); /* Ti */
	P(G, C) = c * (1 - e);                                     /* Tv */
	P(G, G) = (g * (a + g + (c + t) * e) + a * eAG) / (a + g); /* self */
	P(G, T) = t * (1 - e);                                     /* Tv */

	P(T, A) = a * (1 - e);                                     /* Tv */
	P(T, C) = (c * (c + t + (a + g) * e) - c * eCT) / (c + t); /* Ti */
	P(T, G) = g * (1 - e);                                     /* Tv */
	P(T, T) = (t * (c + t + (a + g) * e) + c * eCT) / (c + t); /* self */

	return P;
}

} /* namespace EGriceLab */

#endif /* SRC_HKY85_H_ */
