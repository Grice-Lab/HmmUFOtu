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

}; /* end class StringUtils */
} /* namespace EGriceLab */
#endif /* STRINGUTILS_H_ */
