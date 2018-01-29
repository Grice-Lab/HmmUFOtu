/*
 * bHMM_IO_test.cpp
 *
 *  Created on: Mar 25, 2016
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
		cerr << "Usage:  " << argv[0] << " HMM-INFILE HMM-OUTFILE" << endl;
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

	BandedHMMP7 hmm;
	in >> hmm;
	if(in.good())
		cerr << "Hmm read" << endl;
	else {
		cerr << "Unable to read Hmm file: " << argv[1] << endl;
		return EXIT_FAILURE;
	}

	out << hmm;
	if(out.good())
		cerr << "Hmm written" << endl;
	else {
		cerr << "Unable to write to Hmm file: " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	return 0;
}
