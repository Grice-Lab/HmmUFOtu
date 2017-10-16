/*
 * Stats.h
 *
 *  Created on: Jul 24, 2015
 *      Author: zhengqi
 *      This header includes many basic statistical functions for HmmUFOtu project
 */

#ifndef STATS_H_
#define STATS_H_

#include <cstddef>
#include <map>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <Eigen/Dense>

namespace EGriceLab {

namespace Math {

using std::map;
using std::vector;
using Eigen::VectorXd;
using Eigen::Vector4d;

/**
 * A template method to found the associated key of the maximum value in a std::map
 * The mapped_type of the map must support strict less (operator<)
 * @param freq  a frequency map
 * @return the key whose associated value is maximum
 */
template <typename K, typename V>
K which_max(map<K, V> freq) {
	assert(!freq.empty());
	K maxKey = freq.begin()->first;
	V maxVal = freq.begin()->second;
	for(typename map<K, V>::const_iterator it = freq.begin(); it != freq.end(); ++it)
		if(maxVal < it->second) {
			maxKey = it->first;
			maxVal = it->second;
		}
	return maxKey;
}

/**
 * A template method to found the maximum index in a std::vector
 * The template type T must support strict less (operator<)
 * @param count  a vector of count or frequency
 * @return the index whose value is maximum
 */
template <typename T>
typename vector<T>::size_type which_max(vector<T> count) {
	assert(!count.empty());
	typename vector<T>::size_type idx = 0;
	T max = count[0];
	for(typename vector<T>::size_type i = 1; i != count.size(); ++i)
		if(max < count[i]) {
			idx = i;
			max = count[i];
		}
	return idx;
}

/**
 * A template method to found the maximum index in an array
 * The template type T must support strict less (operator<)
 * @param arr  array to search
 * @param n  array size, must be non-zero
 * @return the index whose value is maximum
 */
template <typename T>
size_t which_max(const T* arr, size_t n) {
	assert(n > 0);
	size_t idx = 0;
	T max = arr[0];
	for(size_t i = 1; i < n; ++i)
		if(max < arr[i]) {
			idx = i;
			max = arr[i];
		}
	return idx;
}

/**
 * A template method to found the maximum value in a std::map
 * The mapped_type of the map must support strict less (operator<)
 * @param freq  a frequency map
 * @return the maximum value
 */
template <typename K, typename V>
V max(map<K, V> freq) {
	assert(!freq.empty());
	V max = freq.begin()->second;
	for(typename map<K, V>::const_iterator it = freq.begin(); it != freq.end(); ++it)
		if(it->second > max) {
			max = it->second;
		}
	return max;
}

/**
 * A template method to found the maximum value in an array
 * The mapped_type of the map must support strict less (operator<)
 * @param arr  array to search
 * @param n  array size, must be non-zero
 * @return the maximum value
 */
template <typename T>
T max(const T* arr, size_t n) {
	assert(n > 0);
	const T* ptr = arr;
	T max = ptr++[0];
	while(ptr != arr + n) {
		if(max < *ptr)
			max = *ptr;
		ptr++;
	}
	return max;
}

/**
 * A template method to check whether a value x is in a given vector
 * The mapped_type of the map must support comparison (operator=)
 * @param x  value to be checked
 * @param vec  vector to be checked
 * @return true if any element in vec equals to x
 */
template <typename T>
bool is_element(T x, vector<T> vec) {
	for(typename vector<T>::const_iterator it = vec.begin(); it != vec.end(); ++it)
		if(*it == x)
			return true;
	return false;
}

/**
 * A template method to count how many times a value x is in a given array
 * The mapped_type of the map must support comparison (operator=)
 * @param x  value to be checked
 * @param array  array to be checked
 * @param n  array size
 * @return number of times x is in array
 */
template <typename T>
size_t count_element(T x, const T* arr, size_t n) {
	size_t count = 0;
	for(const T* ptr = arr; ptr != arr + n; ++ptr)
		if(*ptr == x)
			count++;
	return count;
}

/**
 * A template method to count how many times a value x is in not a given array
 * The mapped_type of the map must support comparison (operator!=)
 * @param x  value to be checked
 * @param array  array to be checked
 * @param n  array size
 * @return number of times x is not in array
 */
template <typename T>
size_t count_not_element(T x, const T* arr, size_t n) {
	size_t count = 0;
	for(const T* ptr = arr; ptr != arr + n; ++ptr)
		if(*ptr != x)
			count++;
	return count;
}

/**
 * calculate bit-per-element for an given alphabet
 * @param n  maximum number to encode
 * @return bits required to encode this numbers upto this value, or -1 if size is zero or negative
 */
inline int bpe(int n) {
	if(n <= 0)
		return -1;
	int shift = 0;
	for(int x = 1; x < n; x <<= 1)
		shift++;
	return shift;
}

/**
 * A template method to calculate the sum of an array
 * The mapped_type of the map must support operator +=
 * @param arr  array
 * @param n  array size
 * @return the sum of the array
 */
template <typename T>
T sum(const T* arr, size_t n) {
	T sum = 0;
	for(const T* ptr = arr; ptr != arr + n; ++ptr)
		sum += *ptr;
	return sum;
}

/**
 * A template method to calculate the weighted sum of an array
 * The mapped_type of the map must support operator +=
 * @param arr  array
 * @param w  weight
 * @param n  array size
 * @return the sum of the array
 */
template <typename T>
double sum(const T* arr, const double* w, size_t n) {
	double sum = 0;
	for(size_t i = 0; i != n; ++i)
		sum += arr[i] * w[i];
	return sum;
}

/**
 * normalize a given double array
 */
inline void normalize(double* arr, size_t n, double C = 1.0) {
	if(n == 0)
		return;
	double s = sum(arr, n);
	for(double* ptr = arr; ptr != arr + n; ++ptr)
		*ptr /= s * C;
}

/**
 * add two log-likelihood value without underflowing by scaling
 */
inline double add_scaled(double logA, double logB) {
	double scale = std::max(logA, logB); /* always scaling */
	return ::log(::exp(logA - scale) + ::exp(logB - scale)) + scale;
}

/** calculate Q value from p-value */
inline double p2q(double p, double b = 10) {
	return -b * ::log(p) / ::log(b);
}

/** calculate p value from q-value */
inline double q2p(double q, double b = 10) {
	return ::exp(- q / b * ::log(b));
}

/** calculate Eucledian distance between two vectors */
inline double euclideanDist(const VectorXd& p, const VectorXd& q) {
	assert(p.rows() == q.rows());
	return (p - q).norm();
}

/** calculate Eucledian distance between two vectors */
inline double euclideanDist(const Vector4d& p, const Vector4d& q) {
	return (p - q).norm();
}

/** calculate Bhattacharyya distance between two vectors */
inline double bhattacharyyaDist(const VectorXd& p, const VectorXd& q) {
	assert(p.rows() == q.rows());
	return - ::log(p.cwiseProduct(q).cwiseSqrt().sum());
}

/** calculate the Kullback–Leibler KL divergence */
inline double KLDivergence(const Vector4d& p, const Vector4d& q) {
	return (p.array() * (p.array() / q.array()).log()).sum();
}

/** calculate the Kullback–Leibler KL divergence */
inline double KLDivergence(const VectorXd& p, const VectorXd& q) {
	return (p.array() * (p.array() / q.array()).log()).sum();
}

} /* namespace Math */

} /* namespace EGriceLab */

#endif /* STATS_H_ */
