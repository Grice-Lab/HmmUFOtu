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
#include "StringUtils.h"

namespace EGriceLab {

string StringUtils::remove_dup_chars(const string& str) {
	string newStr;
	for(string::const_iterator it = str.begin(); it != str.end(); ++it)
		if(newStr.find(*it) == string::npos) // not exist
			newStr.push_back(*it);
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
	if(suffix.empty())
		return StringUtils::basename(path);
	if(suffix.front() != '.')
		suffix.insert(suffix.begin(), '.');

	/* trim directory path */
	path.erase(0, path.find_last_of('/') + 1); /* erase prefix, could be empty (0 length) */
	if(path.length() > suffix.length() && path.substr(path.length() - suffix.length()) == suffix) /* suffix exists */
		path.erase(path.length() - suffix.length());
	return path;
}

string StringUtils::basename(string path) {
	path.erase(0, path.find_last_of('/') + 1); /* erase prefix, could be empty (0 length) */
	string::size_type sufStart = path.find_last_of('.');
	if(sufStart != string::npos)
		path.erase(sufStart);
	return path;
}

string StringUtils::stripQuotes(const string& str, const string& quotes) {
	string newStr;
	newStr.reserve(str.length());
	for(string::const_iterator it = str.begin(); it != str.end(); ++it) {
		if((it == str.begin() || it == str.end() - 1) && /* leading or tailing character */
				quotes.find(*it) != string::npos) /* is a quote character */
			continue;
		newStr.push_back(*it);
	}
	return newStr;
}

bool StringUtils::containsWhiteSpace(const string& str) {
	for(string::const_iterator it = str.begin(); it != str.end(); ++it)
		if(::isspace(*it))
			return true;
	return false;
}

bool StringUtils::containsAny(const string& str, const string& query) {
	for(string::const_iterator it = query.begin(); it != query.end(); ++it)
		if(str.find(*it) != string::npos)
			return true;
	return false;
}

string& StringUtils::removeAll(string& str, const string& pattern) {
	string::size_type n = pattern.length();
	for(string::size_type i = str.find(pattern); i!= string::npos; i = str.find(pattern))
		str.erase(i, n);
	return str;
}

string& StringUtils::removeEnd(string& str, const string& suffix) {
	if(str.rfind(suffix) == str.length() - suffix.length())
		str.erase(str.end() - suffix.length(), str.end());
	return str;
}

string::size_type StringUtils::common(const string& str1, const string& str2) {
	string::size_type N = 0;
	string::size_type count1[CHAR_MAX + 1] = { }; /* zero initialization */
	string::size_type count2[CHAR_MAX + 1] = { }; /* zero initialization */

	for(string::const_iterator it = str1.begin(); it != str1.end(); ++it)
		count1[*it]++;
	for(string::const_iterator it = str2.begin(); it != str2.end(); ++it)
		count2[*it]++;
	for(int i = 0; i <= CHAR_MAX; ++i)
		if(count1[i] && count2[i])
			N++;
	return N;
}

size_t StringUtils::common(const char* str1, const char* str2) {
	size_t N = 0;
	size_t count1[CHAR_MAX + 1] = { }; /* zero initialization */
	size_t count2[CHAR_MAX + 1] = { }; /* zero initialization */

	for(; *str1; ++str1)
		count1[*str1]++;
	for(; *str2; ++str2)
		count2[*str2]++;
	for(int i = 0; i <= CHAR_MAX; ++i)
		if(count1[i] && count2[i])
			N++;
	return N;
}

} /* namespace EGriceLab */

