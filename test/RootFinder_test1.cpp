/*
 * RootFinder_test1.cpp
 *  Unit test for RootFinder class
 *  Created on: Oct 24, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "math/RootFinder.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::Math;

/**
 * A concrete subclass of RootFinder::R2RFunc to solve quadratic function as
 * f(x) = ax^2 + bx + c = 0
 */
struct QuadFunc : RootFinder::R2RFunc {
	QuadFunc(double a, double b, double c) : a(a), b(b), c(c) { }

	virtual ~QuadFunc() { }

	virtual double operator()(double x) {
		return a * x * x + b * x + c;
	}

	double a, b, c;
};

int main(int argc, char *argv[]) {
	if(argc != 6) {
		cerr << "Usage:  " << argv[0] << " a  b  c xl xr" << endl;
		return -1;
	}
	int a = ::atoi(argv[1]);
	int b = ::atoi(argv[2]);
	int c = ::atoi(argv[3]);

	double xl = ::atof(argv[4]);
	double xr = ::atof(argv[5]);

	QuadFunc qf(a, b, c);
	cerr << "Quadratic function " << a << "x^2 + " << b << "x + " << c << " constructed" << endl;

	RootFinder rf(qf, xl, xr);
	cerr << "Finding root in [" << xl << ", " << xr << "]" << endl;

	double x = rf.rootBisection();
	cerr << "Root found at: " << x << " where f(x) = " << qf(x) << endl;
}



