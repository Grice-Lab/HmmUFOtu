/*
 * alphabet_test1.cpp
 *
 *  Created on: Oct 26, 2015
 *      Author: zhengqi
 */


#include <iostream>
#include <cstdlib>
#include <ctime>
#include "IUPACNucl.h"

using namespace std;
using namespace EGriceLab;

int main() {
	::srand(time(NULL));
	cout << "Initiating a test IUPACNucl object" << endl;
	IUPACNucl nucl1;
	for(string::const_iterator it = nucl1.getSymbol().begin(); it != nucl1.getSymbol().end(); ++it)
		cout << *it << " -> " << static_cast<int> (nucl1.encode(*it)) << endl;

	for(string::const_iterator it = nucl1.getSynonymous().begin(); it != nucl1.getSynonymous().end(); ++it)
		cout << *it << " -> " << static_cast<int> (nucl1.encode(*it)) << endl;

	cout << "Copying the test IUPACNucl pointer" << endl;
	const DegenAlphabet* nucl2 = new IUPACNucl(nucl1);

	for(string::const_iterator it = nucl2->getSymbol().begin(); it != nucl2->getSymbol().end(); ++it)
		cout << *it << " -> " << static_cast<int> (nucl2->encode(*it)) << endl;
	for(string::const_iterator it = nucl2->getSynonymous().begin(); it != nucl2->getSynonymous().end(); ++it)
		cout << *it << " -> " << static_cast<int> (nucl2->encode(*it)) << endl;

	cout << "nucl1 == nucl2:" << (nucl1 == *nucl2) << endl;

	return 0;
}
