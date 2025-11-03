/*
 * StringUtil_test.cpp
 *
 *  Created on: May 21, 2021
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstdlib>
#include "StringUtils.h"
using namespace std;
using namespace EGriceLab;

int main() {
	const string str1 = " \tHello World\t ";
	const string str2 = "\"'Hello World'\"";
	const string strTrimmed = "Hello World";

	cout << "str1:" << endl << str1 << endl;
	cout << "str2:" << endl << str2 << endl;

	if(StringUtils::trim(str1) != strTrimmed) {
		cerr << "Unmatched result:" << endl << StringUtils::trim(str1) << endl << strTrimmed << endl;
		return EXIT_FAILURE;
	}

	if(StringUtils::stripQuotes(str2) != strTrimmed) {
		cerr << "Unmatched result:" << endl << StringUtils::trim(str2) << endl << strTrimmed << endl;
		return EXIT_FAILURE;
	}
	return 0;
}
