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
		b(i) = quantile(gammaDist, i / static_cast<double>(K));
	b(K) = inf;
}

void DiscreteGammaModel::setRates() {
	for(int i = 0; i < K; ++i) {
		double lbd = b(i);
		double ubd = b(i+1);
		r(i) = ubd != inf ?
				gamma_p(alpha + 1, ubd * alpha) - gamma_p(alpha + 1, lbd * alpha) :
				1 - gamma_p(alpha + 1, lbd * alpha);
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

	map.segment(0, K) = r; /* copy r into buf */
	out.write((const char*) buf, sizeof(double) * K);
	delete[] buf;

	return out;
}

double DiscreteGammaModel::estimateShape(const VectorXi& X) {
	if(X.rows() < 2)
		return EGriceLab::inf; // cannot estimate alpha, use inf
	double m = X.mean();
	double s = (X.array() - m).matrix().squaredNorm() / (X.rows() - 1);
	return m * m / (s - m);
}

} /* namespace EGriceLab */
