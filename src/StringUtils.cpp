/*
 * StringUtils.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: zhengqi
 */

#include <algorithm>
#include <cctype>
#include <iostream>
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

string StringUtils::toLower(const string& str) {
	string newStr; // make a new copy
	newStr.resize(str.length());
	transform(str.begin(), str.end(), newStr.begin(), ::tolower);
	return newStr;
}

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
	/* trim directory path*/
	path.erase(0, path.find_last_of('/') + 1); /* erase prefix, could be empty (0 length) */
	if(!suffix.empty()) { /* suffix specified */
		if(suffix[0] != '.')
			suffix = '.' + suffix;
		if(path.length() > suffix.length() && path.substr(path.length() - suffix.length()) == suffix) /* suffix exists */
			path.erase(path.length() - suffix.length());
	}
	else { /* trim any suffix */
		size_t sufStart = path.find_last_of('.');
		if(sufStart != string::npos)
			path.erase(sufStart);
	}
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

} /* namespace EGriceLab */
