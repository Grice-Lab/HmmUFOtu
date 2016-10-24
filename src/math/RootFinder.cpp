/*
 * RootFinder.cpp
 *
 *  Created on: Oct 21, 2016
 *      Author: zhengqi
 */

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include "RootFinder.h"

namespace EGriceLab {
namespace Math {

//const double RootFinder::DEFAULT_ABS_EPS = FLT_EPSILON;
//const double RootFinder::DEFAULT_REL_EPS = FLT_EPSILON;
const double RootFinder::DEFAULT_ABS_EPS = 1e-10;
const double RootFinder::DEFAULT_REL_EPS = 1e-10;
const double RootFinder::DEFAULT_RES_EPS = 0;

double RootFinder::rootBisection() {
	double x, xmag, fx;
	int iter;

	for(iter = 0; maxIter == 0 || iter < maxIter; ++iter) {
		/* Bisect and evaluate the function */
		x = (xl + xr) / 2;
		fx = f(x);
		if(fx == 0) /* an exact root, lucky */
			break;

		/* test for convergence */
		double fxl = f(xl);
		double fxr = f(xr);
		xmag = (xl < 0 && xr > 0) ? 0 : x;

		if(xr - xl < absEps + relEps * xmag || ::fabs(fx) < resEps) /* an approximate root */
			break;

		/* narrow the bracket */
		if(fxl > 0) {
			if(fx > 0) {
				xl = x;
				fxl = fx;
			}
			else {
				xr = x;
				fxr = fx;
			}
		}
		else {
			if(fx < 0) {
				xl = x;
				fxl = fx;
			}
			else {
				xr = x;
				fxr = fx;
			}
		}
	}

	if(maxIter > 0 && iter >= maxIter) {
		std::cerr << "RootFinder unable to converge after " << maxIter << " iteration" << std::endl;
		std::abort();
	}
	return x;
}

} /* namespace Math */
} /* namespace EGriceLab */

