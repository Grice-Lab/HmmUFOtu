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
 * StringUtils.h
 *
 *  Created on: Jul 22, 2015
 *      Author: zhengqi
 *  This header contains non-class utility functions to handle common string manipulations
 */

#ifndef STRINGUTILS_H_
#define STRINGUTILS_H_
#include <string>
#include <iostream>

namespace EGriceLab {
using std::string;
using std::basic_string;
using std::istream;
using std::ostream;

class StringUtils {
public:
	/**
	 * remove duplicate characters in a string
	 * @param str  target string
	 * @return a copy with duplicated characters removed
	 */
	static string remove_dup_chars(const string& str);

	/**
	 * make a string to all uppercase
	 * @param str  original string
	 * @return a new copy with all upper-case as the original string
	 */
	static string toUpper(const string& str);

	/**
	 * make a string to all uppercase
	 * @param str  original string
	 * @return the modified string with characters in all upper-case as the original string
	 */
	static string& toUpper(string& str);

	/**
	 * make a string to all lower-case
	 * @param str  original string
	 * @return a new copy with all upper-case as the original string
	 */
	static string toLower (const string& str);

	/**
	 * make a string to all lower-case
	 * @param str  original string
	 * @return the modified string with characters in all upper-case as the original string
	 */
	static string& toLower(string& str);

	/**
	 * Test whether a string ends with a given suffix
	 * @param str  target string
	 * @param suffix  suffix to check
	 * @return  true only if str ends with suffix
	 */
	static bool endsWith(const string& str, const string& suffix);

	/**
	 * Test whether a string starts with a given prefix
	 * @param str  target string
	 * @param prefix  suffix to check
	 * @return  true only if str starts with suffix
	 */
	static bool startsWith(const string& str, const string& prefix);

	/**
	 * Get basename of a path/filename, trim all leading path and optionally tailing suffix, if any
	 * @param @path  pathname
	 * @param @suffix  suffix of filename
	 * @return a new string with directory path and suffix trimmed
	 */
	static string basename(string path, string suffix);

	/**
	 * Get basename of a path/filename, trim all leading path tailing suffix, if any
	 * @param @path  pathname
	 * @return a new string with directory path and suffix trimmed
	 */
	static string basename(string path);

	/**
	 * Remove leading and tailing quotes from a given string
	 * @param str  string input
	 * @param quotes  quoting characters
	 * @return a new string with all quotes in "quote" removed
	 */
	static string stripQuotes(const string& str, const string& quotes = "\"'");

	/**
	 * check whether a string contains any white space characters
	 * @param str  input string
	 * @return  true if it has any white space character (' ', '\t', '\n', '\r', '\v')
	 */
	static bool containsWhiteSpace(const string& str);

	/**
	 * check whether a string contains any character in another string
	 * @param str  input string
	 * @param query  query string
	 * @return  true if str contains any character of query
	 */
	static bool containsAny(const string& str, const string& query);

	/**
	 * Remove all occurrences of pattern in str
	 * @param str  input string
	 * @param pattern  pattern to remove
	 * @return  the modified input
	 */
	static string& removeAll(string& str, const string& pattern);

	/**
	 * Remove the given tail in a string, if exists
	 * @param str  input string
	 * @param suffix  siffix to remove
	 * @return  the modified input
	 */
	static string& removeEnd(string& str, const string& pattern);

	/**
	 * load data from a binary input to given basic_string, override any old data
	 * @param dest  destination
	 * @param in  input
	 * @param number basic_string to load
	 * @return  whether loading was successful
	 */
	template<typename T>
	static bool loadString(basic_string<T>& dest, istream& in, size_t length) {
		T* buf = new T[length]; /* construct a temporary buffer */
		in.read((char*) buf, length * sizeof(T));
		dest.assign(buf, length);
		delete[] buf;

		return in.good();
	}

	/**
	 * save a basic_string to a binary output, upto length of source will be saved
	 * @param src  source
	 * @param out  output
	 * @param length  number of source to save
	 * @return  whether saving was successful
	 */
	template<typename T>
	static bool saveString(const basic_string<T>& src, ostream& out, size_t length) {
		out.write((const char*) src.c_str(), length * sizeof(T));
		return out.good();
	}

	template<typename T>
	static bool saveString(const basic_string<T>& src, ostream& out) {
		return saveString(src, out, src.length());
	}

	static string::size_type common(const string& str1, const string& str2);

	static size_t common(const char* str1, const char* str2);

}; /* end class StringUtils */

} /* namespace EGriceLab */
#endif /* STRINGUTILS_H_ */
