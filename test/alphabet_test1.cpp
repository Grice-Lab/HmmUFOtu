/*
 * alphabet_test1.cpp
 *
 *  Created on: Oct 26, 2015
 *      Author: zhengqi
 */



#include <iostream>
#include "Alphabet.h"

using namespace std;
using namespace EGriceLab;

int main() {

	cout << "Initiating a test DNA alphabet" << endl;
	Alphabet abc("DNA", "ATCG");

	for(string::const_iterator it = abc.getSymbol().begin(); it != abc.getSymbol().end(); ++it)
		cout << *it << " -> " << static_cast<int> (abc.encode(*it)) << endl;
	return 0;
}
