/*
 * MultinomialDistribution_test.cpp
 *
 *  Created on: May 14, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include "MultinomialDistribution.h"

using namespace std;
using namespace Eigen;
using namespace EGriceLab::Math;

int main() {
	VectorXd p = VectorXd::Ones(4) / 4.0;
	cout << "p: " << p.transpose() << endl;
	multinomial dist(p);

	VectorXd x(4);

	x << 25.0, 25.0, 25.0, 25.0;
	cout << "x: " << x.transpose() << endl;
	cout << "pdf(x, p): " << pdf(dist, x) << endl;

	x << 20.0, 30.0, 15.0, 35.0;
	cout << "x: " << x.transpose() << endl;
	cout << "pdf(x, p): " << pdf(dist, x) << endl;

	x << 10.0, 40.0, 5.0, 45.0;
	cout << "x: " << x.transpose() << endl;
	cout << "pdf(x, p): " << pdf(dist, x) << endl;

	VectorXf q = VectorXf::Ones(4) / 4.0;
	cout << "q: " << q.transpose() << endl;
	MultinomialDistribution<float> dist2(q);
	VectorXf y(4);

	y << 25.0, 25.0, 25.0, 25.0;
	cout << "y: " << y.transpose() << endl;
	cout << "pdf(y, q): " << pdf(dist2, y) << endl;

	y << 20.0, 30.0, 15.0, 35.0;
	cout << "y: " << y.transpose() << endl;
	cout << "pdf(y, q): " << pdf(dist2, y) << endl;

	y << 10.0, 40.0, 5.0, 45.0;
	cout << "y: " << y.transpose() << endl;
	cout << "pdf(y, q): " << pdf(dist2, y) << endl;
}
