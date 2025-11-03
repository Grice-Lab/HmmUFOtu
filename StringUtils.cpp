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
 * StringUtils.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: zhengqi
 */

#include <algorithm>
#include <cctype>
#include <iostream>
#include <climits>
#include <array>
#include <regex>
#include "StringUtils.h"

namespace EGriceLab {

string StringUtils::remove_dup_chars(const string& str) {
	string newStr;
	for(string::value_type c : str)
		if(newStr.find(c) == string::npos) // not exist
			newStr.push_back(c);
	return newStr;
}

string StringUtils::toUpper(const string& str) {
	string newStr; // make a new copy
	newStr.resize(str.length());
	transform(str.begin(), str.end(), newStr.begin(), ::toupper);
	return newStr;
}

string& StringUtils::toUpper(string& str) {
	transform(str.begin(), str.end(), str.begin(), ::toupper);
	return str;
}

/**
 * make a copy of the input string in all lower cases
 * @param str  input string
 * @return  a copy with in all lower cases
 */
string StringUtils::toLower(const string& str) {
	string newStr; // make a new copy
	newStr.resize(str.length());
	transform(str.begin(), str.end(), newStr.begin(), ::tolower);
	return newStr;
}

/**
 * make the input string into all lower cases
 * @param str  input string
 * @return  the modified string
 */
string& StringUtils::toLower(string& str) {
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	return str;
}

bool StringUtils::endsWith(const string& str, const string& suffix) {
	if(str.length() < suffix.length())
		return false;
	return str.substr(str.length() - suffix.length()) == suffix;
}

bool StringUtils::startsWith(const string& str, const string& prefix) {
	if(str.length() < prefix.length())
		return false;
	return str.substr(0, prefix.length()) == prefix;
}

string StringUtils::basename(string path, string suffix) {
	/* trim directory path */
	path.erase(0, path.find_last_of('/') + 1); /* erase prefix, could be empty (0 length) */
	/* trim optional suffix */
	if(!suffix.empty()) {
		if(suffix[0] != '.')
			suffix.insert(suffix.begin(), '.');
		if(path.length() > suffix.length() && path.substr(path.length() - suffix.length()) == suffix) /* suffix exists */
			path.erase(path.length() - suffix.length());
	}
	return path;
}

string StringUtils::trim(const string& str, const string& ws) {
	string newStr;
	newStr.reserve(str.length());
	/* construct re */
	std::regex re("(?:^[" + ws + "]+)|(?:[" + ws + "]+$)");
	/* strip heading/tailing quotes */
	newStr = std::regex_replace(str, re, "");
	return newStr;
}

string StringUtils::stripQuotes(const string& str, const string& quotes) {
	string newStr;
	newStr.reserve(str.length());
	/* construct regex */
	std::regex re("(?:^[" + quotes + "]+)|(?:[" + quotes + "]+$)");
	/* strip heading/tailing quotes */
	newStr = std::regex_replace(str, re, "");
	return newStr;
}

bool StringUtils::containsWhiteSpace(const string& str) {
	return std::any_of(str.begin(), str.end(), ::isspace);
}

bool StringUtils::containsAny(const string& str, const string& query) {
	return std::any_of(query.begin(), query.end(),
			[=] (string::value_type c) { return str.find(c) != string::npos; });
}

string& StringUtils::removeAll(string& str, const string& pattern) {
	string::size_type n = pattern.length();
	for(string::size_type i = str.find(pattern); i!= string::npos; i = str.find(pattern))
		str.erase(i, n);
	return str;
}

string StringUtils::removeAll(const string& str, const string& pattern) {
	string strN = str;
	string::size_type n = pattern.length();
	for(string::size_type i = strN.find(pattern); i!= string::npos; i = strN.find(pattern))
		strN.erase(i, n);
	return strN;
}

string& StringUtils::removeEnd(string& str, const string& suffix) {
	if(str.rfind(suffix) == str.length() - suffix.length())
		str.erase(str.end() - suffix.length(), str.end());
	return str;
}

string StringUtils::removeEnd(const string& str, const string& suffix) {
	string strN = str;
	if(strN.rfind(suffix) == strN.length() - suffix.length())
		strN.erase(strN.end() - suffix.length(), strN.end());
	return strN;
}

string StringUtils::common(string str1, string str2) {
	/* sort input strings */
	std::sort(str1.begin(), str1.end());
	std::sort(str2.begin(), str2.end());
	string::size_type n = std::min(str1.length(), str2.length());
	string comm(n, '\0'); // construct a common string with enough space
	string::iterator result = std::set_intersection(str1.begin(), str1.end(),
			str2.begin(), str2.end(), comm.begin());
	comm.resize(result - comm.begin()); // resize to fit
	return comm;
}

} /* namespace EGriceLab */

