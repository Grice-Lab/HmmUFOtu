/*
 * RootFinder.h
 *  Finding roots of a functor that map in R -> R (real) space
 *  Created on: Oct 21, 2016
 *      Author: zhengqi
 */

#ifndef SRC_MATH_ROOTFINDER_H_
#define SRC_MATH_ROOTFINDER_H_

#include <stdexcept>

namespace EGriceLab {
namespace Math {

using std::invalid_argument;

class RootFinder {
public:
	/* enclosing types and enums */
	/**
	 * An abstracted functor operator()(x) T->R functor that can be used to find the root by chaning the value of x
	 */
	struct R2RFunc {
		/*
		 * pure virtual (abstract) method
		 * evaluate the functor at x
		 * @param x  functor parameter
		 * @return  result evaluated at x
		 */
		virtual double operator()(double x) = 0;

		/**
		 * virtual destructor, do nothing
		 */
		virtual ~R2RFunc() { }
	};

	/* constructors and destructor */
	/**
	 * construct a RootFinder in given domain [xl, xr]
	 *
	 */
	RootFinder(R2RFunc& f, double xl, double xr) :
		f(f), xl(xl), xr(xr),
		absEps(DEFAULT_ABS_EPS), relEps(DEFAULT_REL_EPS),
		resEps(DEFAULT_RES_EPS), maxIter(MAX_ITER) {
		if(setDomain(xl, xr) >= 0)
			throw invalid_argument("xl, xr do not bracket a root");
	}

	/* disable copy and assignment constructor */
private:
	RootFinder(const RootFinder& other);
	RootFinder& operator=(const RootFinder& other);

public:
	/* member methods */
	/**
	 * set the root search domain
	 * @param xl  lower search bound
	 * @param xr  upper search bound
	 * return  f(xl) * f(xr)
	 */
	double setDomain(double xl, double xr) {
		this->xl = xl;
		this->xr = xr;
		return f(xl) * f(xr);
	}

	/**
	 * Set absolute epsilon
	 */
	void setAbsEps(double absEps) {
		this->absEps = absEps;
	}

	/**
	 * Set relative epsilon
	 */
	void setRelEps(double relEps) {
		this->relEps = relEps;
	}

	/**
	 * Set residue epsilon
	 */
	void setResEps(double resEps) {
		this->resEps = resEps;
	}

	/**
	 * Set maximum iteration
	 */
	void setMaxIter(int maxIter) {
		this->maxIter = maxIter;
	}

	/**
	 * find one-dimensional root of the functor f in a new domain
	 * @param xl  lower search bound
	 * @param xr  upper search bound
	 * return root x so f(x) == 0
	 */
	double rootBisection(double xl, double xr);

	/**
	 * find one-dimensional root of the functor f using current domain
	 * return root x so f(x) == 0
	 */
	double rootBisection();

private:
	R2RFunc& f;

	double xl;
	double xr;
//	double fl;
//	double fr;
//	double x;
//	double fx;
	double absEps;
	double relEps;
	double resEps;
	int maxIter;

	static const double DEFAULT_ABS_EPS; /* absolute epsilon */
	static const double DEFAULT_REL_EPS; /* relative epsilon */
	static const double DEFAULT_RES_EPS; /* residue epsilon */
	static const int MAX_ITER = 0;       /* maximum iteration */
};


inline double RootFinder::rootBisection(double xl, double xr) {
	if(setDomain(xl, xr) >= 0)
		throw invalid_argument("xl, xr do not bracket the root");
	return rootBisection();
}

} /* namespace Math */
} /* namespace EGriceLab */

#endif /* SRC_MATH_ROOTFINDER_H_ */
