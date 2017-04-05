/*
 * TN93.h
 *  TN93 DNA Substitution Model
 *  Created on: Mar 8, 2017
 *      Author: zhengqi
 */

#ifndef SRC_TN93_H_
#define SRC_TN93_H_

#include <cmath>
#include "DNASubModel.h"

namespace EGriceLab {

class TN93: public DNASubModel {
public:
	/* Constructors */

	/** default constructor */
	TN93() : kr(1), ky(1), pi(Vector4d::Constant(1.0/4))
	{
		setBeta();
	}

	/* destructor, do nothing */
	virtual ~TN93() { }

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
	 * The formular is discribed in the original TN93 article
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
	virtual TN93* clone() const {
		return new TN93(*this);
	}

private:
	/** set beta by kappa and pi */
	void setBeta() {
		beta = 1 / ( 2 * (pi(A) * pi(C) + pi(A) * pi(T) + pi(C) * pi(G) + pi(G) * pi(T)
				       + kr * (pi(A) * pi(G)) + ky * (pi(C) * pi(T))));
	}

	static const string name;

	Vector4d pi; /* base frequency */
	double kr; // Ti/Tv ratio of purines
	double ky; // Ti/Tv ratio of pyrimidines
	double beta; // sequence diversity as 1 / 2(AC + AT + CG + GT + kr * AG + ky * CT))
};

inline Matrix4d TN93::Pr(double v) const {
	Matrix4d P;
	double a = pi(A);
	double c = pi(C);
	double g = pi(G);
	double t = pi(T);
	double e = ::exp(-beta * v);
	double eR = ::exp(-(1 + (a + g) * (kr - 1)) * beta * v);
	double eY = ::exp(-(1 + (c + t) * (ky - 1)) * beta * v);

	P(A, A) = (a * (a + g + (c + t) * e) + g * eR) / (a + g); /* self */
	P(A, C) = c * (1 - e);                                    /* Tv */
	P(A, G) = (g * (a + g + (c + t) * e) - g * eR) / (a + g); /* Ti */
	P(A, T) = t * (1 - e);                                    /* Tv */

	P(C, A) = a * (1 - e);                                    /* Tv */
	P(C, C) = (c * (c + t + (a + g) * e) + t * eY) / (c + t); /* self */
	P(C, G) = g * (1 - e);                                    /* Tv */
	P(C, T) = (t * (c + t + (a + g) * e) - t * eY) / (c + t); /* Ti */

	P(G, A) = (a * (a + g + (c + t) * e) - a * eR) / (a + g); /* Ti */
	P(G, C) = c * (1 - e);                                    /* Tv */
	P(G, G) = (g * (a + g + (c + t) * e) + a * eR) / (a + g); /* self */
	P(G, T) = t * (1 - e);                                    /* Tv */

	P(T, A) = a * (1 - e);                                    /* Tv */
	P(T, C) = (c * (c + t + (a + g) * e) - c * eY) / (c + t); /* Ti */
	P(T, G) = g * (1 - e);                                    /* Tv */
	P(T, T) = (t * (c + t + (a + g) * e) + c * eY) / (c + t); /* self */

	return P;
}

inline double TN93::subDist(const Matrix4d& D, double N) const {
	if(N == 0)
		return 0;
	double a = pi(A);
	double c = pi(C);
	double g = pi(G);
	double t = pi(T);
	double r = a + g;
	double y = c + t;
	double pr = (D(A,G) + D(G,A)) / N; /* observed Tr diff */
	double py = (D(C,T) + D(T,C)) / N; /* observed Ty diff */
	double q = (D(A,C) + D(A,T) + D(C,A) + D(C,G) + D(G,C) + D(G,T) + D(T,A) + D(T,G)) / N; /* observed Tv diff */

	return - 2 * a * g / r * ::log(1 - r / (2 * a * g) * pr - 1 / (2 * r) * q)
		   - 2 * g * c / y * ::log(1 - y / (2 * t * c) * py - 1 / (2 * y) * q)
	       - 2 * (r * y - a * g * y / r - t * c * r / y) * ::log(1 - 1 / (2 * r * y) * q);
}

} /* namespace EGriceLab */

#endif /* SRC_TN93_H_ */
