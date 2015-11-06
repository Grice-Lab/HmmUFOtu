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
#include <cassert>

namespace EGriceLab {
using std::map;
using std::vector;

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
	T max = arr++[0];
	while(arr != arr + n) {
		if(max < *arr)
			max = *arr;
		arr++;
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

} /* namespace EGriceLab */

#endif /* STATS_H_ */
