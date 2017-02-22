/*
 * DiscreteGammaModel.cpp
 *
 *  Created on: Feb 17, 2017
 *      Author: zhengqi
 */

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/unordered_map.hpp>
#include "DiscreteGammaModel.h"
using boost::math::gamma_distribution;
using boost::math::quantile;
using boost::math::gamma_p;
using boost::unordered_map;
using Eigen::Map;

namespace EGriceLab {

void DiscreteGammaModel::setBreaks() {
	gamma_distribution<double> gammaDist(alpha, alpha);
	for(int i = 0; i < K; ++i)
		b(i) = quantile(gammaDist, i / K);
	b(K) = inf;
}

void DiscreteGammaModel::setRates() {
	for(int i = 0; i < K; ++i) {
		double lbd = b(i);
		double ubd = b(i+1);
		r(i) = (gamma_p(ubd * alpha, alpha + 1) - gamma_p(lbd * alpha, alpha + 1)) / K;
	}
}

istream& DiscreteGammaModel::load(istream& in) {
	/* read basic fields */
	in.read((char*) &K, sizeof(int));
	in.read((char*) &alpha, sizeof(double));

	/* read aux fields */
	double *buf = new double[K + 1];
	Map<VectorXd> map(buf, K + 1);
	in.read((char *) buf, sizeof(double) * (K + 1));
	b = map;

	in.read((char*) buf, sizeof(double) * K);
	r = map.segment(0, K); /* ignore the last value */

	return in;
}

ostream& DiscreteGammaModel::save(ostream& out) const {
	/* write basic fields */
	out.write((const char*) &K, sizeof(int));
	out.write((const char*) &alpha, sizeof(double));

	/* write aux fields */
	double *buf = new double[K + 1];
	Map<VectorXd> map(buf, K + 1);
	map = b; /* copy b into buf */
	out.write((const char*) buf, sizeof(double) * (K + 1));

	map = r; /* copy r into buf */
	out.write((const char*) buf, sizeof(double) * K);
	delete[] buf;

	return out;
}

double DiscreteGammaModel::estimateShape(const VectorXi& X) {
	if(X.cols() < 2)
		return EGriceLab::inf; // cannot estimate alpha, use inf
	double m = X.mean();
	double s = (X - m).squaredNorm() / (X.cols() - 1);
	return m * m / (s - m);
}

} /* namespace EGriceLab */
