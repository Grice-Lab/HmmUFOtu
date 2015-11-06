/*
 * alphabet_test1.cpp
 *
 *  Created on: Oct 26, 2015
 *      Author: zhengqi
 */


#include <iostream>
#include <cstdlib>
#include <ctime>
#include "DNA.h"

using namespace std;
using namespace EGriceLab;

int main() {
	::srand(time(NULL));
	cout << "Initiating a test DNA object" << endl;
	DNA dna1;

	for(string::const_iterator it = dna1.getSynonymous().begin(); it != dna1.getSynonymous().end(); ++it)
		cout << *it << " -> " << static_cast<int> (dna1.encode(*it)) << endl;

	cout << "Copying the test DNA pointer" << endl;
	const DNA* dna2 = new DNA(dna1);

	for(string::const_iterator it = dna2->getSynonymous().begin(); it != dna2->getSynonymous().end(); ++it)
		cout << *it << " -> " << static_cast<int> (dna2->encode(*it)) << endl;

	cout << "dna1 == dna2:" << (dna1 == *dna2) << endl;

	return 0;
}
