/*
 * alphabet_test1.cpp
 *
 *  Created on: Oct 26, 2015
 *      Author: zhengqi
 */


#include <iostream>
#include <cstdlib>
#include <ctime>
#include "IUPACAmino.h"

using namespace std;
using namespace EGriceLab;

int main() {
	::srand(time(NULL));
	cout << "Initiating a test IUPACNucl object" << endl;
	IUPACAmino amino1;
	for(string::const_iterator it = amino1.getSymbol().begin(); it != amino1.getSymbol().end(); ++it)
		cout << *it << " -> " << static_cast<int> (amino1.encode(*it)) << endl;

	for(string::const_iterator it = amino1.getSynonymous().begin(); it != amino1.getSynonymous().end(); ++it)
		cout << *it << " -> " << static_cast<int> (amino1.encode(*it)) << endl;

	cout << "Copying the test IUPACNucl pointer" << endl;
	const DegenAlphabet* amino2 = new IUPACAmino(amino1);

	for(string::const_iterator it = amino2->getSymbol().begin(); it != amino2->getSymbol().end(); ++it)
		cout << *it << " -> " << static_cast<int> (amino2->encode(*it)) << endl;
	for(string::const_iterator it = amino2->getSynonymous().begin(); it != amino2->getSynonymous().end(); ++it)
		cout << *it << " -> " << static_cast<int> (amino2->encode(*it)) << endl;

	cout << "amino1 == amino2:" << (amino1 == *amino2) << endl;

	return 0;
}
