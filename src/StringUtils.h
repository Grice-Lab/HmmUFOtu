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

namespace EGriceLab {
using std::string;

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
	static string basename(string path, string suffix = "");

	/**
	 * Remove leading and tailing quotes from a given string
	 * @param str  string input
	 * @param quotes  quoting characters
	 * @return a new string with all quotes in "quote" removed
	 */
	static string stripQuotes(const string& str, const string& quotes = "\"'");

}; /* end class StringUtils */
} /* namespace EGriceLab */
#endif /* STRINGUTILS_H_ */
