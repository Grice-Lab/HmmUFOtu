/*
 * MultinomialDistribution.h
 *
 *  Created on: May 14, 2018
 *      Author: zhengqi
 */

#ifndef MULTINOMIALDISTRIBUTION_H_
#define MULTINOMIALDISTRIBUTION_H_

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <Eigen/Dense>

namespace EGriceLab {
namespace Math {
using Eigen::Matrix;

template <typename RealType = double> class MultinomialDistribution;
typedef MultinomialDistribution<> multinomial;

/** C++ Boost Distribution like class of Multinomial distribution */
template <typename RealType>
class MultinomialDistribution {
public:
	typedef RealType value_type;
	typedef Matrix<RealType, Eigen::Dynamic, 1> VectorXr;

	/* constructors */
	explicit MultinomialDistribution(VectorXr p) : p(p)
	{  }

	/* member methods */
	size_t getK() const {
		return p.size();
	}

	/* member fields */
private:
	VectorXr p;

	/* non-member accessor functions */
	/**
	 * Log-Probability Density (Mass) Function of a Multinomial Distribution
	 * can be expressed as log-gamma function as
	 * lpdf(x, p) = lgamma(sum(x+1)) - sum(lgamma(xi + 1)) + sum(xi * log(pi))
	 */
	friend RealType lpdf(const MultinomialDistribution<RealType>& dist, const VectorXr& x) {
		assert(dist.getK() == x.size());
		/* calculate numerator and denominator */
		RealType num = boost::math::lgamma(x.sum() + 1);
		RealType den = 0;
		for(size_t i = 0; i < x.size(); ++i)
			den += boost::math::lgamma(x(i) + 1);
		/* calculate terms */
		RealType y = (x.cwiseProduct(dist.p.array().log().matrix())).sum();
		return num - den + y;
	}

	/** Probability Density (Mass) Function */
	friend RealType pdf(const MultinomialDistribution<RealType>& dist, const VectorXr& x) {
		return std::exp(lpdf(dist, x));
	}
};

} /* namespace Math */
} /* namespace EGriceLab */

#endif /* MULTINOMIALDISTRIBUTION_H_ */
