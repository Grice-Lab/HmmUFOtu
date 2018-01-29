/*
 * bHmmPrior_IO_test.cpp
 *
 *  Created on: May 31, 2017
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_hmm.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::HmmUFOtu;

int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " PRIOR-INFILE PRIOR-OUTFILE" << endl;
		return 0;
	}

	ifstream in(argv[1]);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	ofstream out(argv[2]);
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	BandedHMMP7Prior hmmPrior;
	in >> hmmPrior;
	if(in.bad()) {
		cerr << "Unable to read HmmPrior file: " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	else
		cerr << "HmmPrior read" << endl;

	out << hmmPrior;
	if(out.good())
		cerr << "HmmPrior written" << endl;
	else {
		cerr << "Unable to write to HmmPrior file: " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	return 0;
}
