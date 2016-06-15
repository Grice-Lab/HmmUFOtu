/*
 * bHMM_IO_test1.cpp
 *
 *  Created on: Mar 25, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include "HmmUFOtu.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " HMM-FILE OUTFILE" << endl;
		return 0;
	}

	ifstream in(argv[1]);
	ofstream out(argv[2]);

	BandedHMMP7 hmm;
	in >> hmm;

	out << hmm;

	in.close();
	out.close();
	return 0;
}
