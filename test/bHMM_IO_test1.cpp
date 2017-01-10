/*
 * bHMM_IO_test1.cpp
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

int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " HMM-FILE OUTFILE" << endl;
		return 0;
	}

	ifstream in(argv[1]);
	if(!in.is_open())
		cerr << "Unable to open " << argv[1] << endl;
	ofstream out(argv[2]);
	if(!out.is_open())
		cerr << "Unable to write to " << argv[2] << endl;

	BandedHMMP7 hmm;
	cerr << "hmm constructed" << endl;
	in >> hmm;
	if(in.good())
		cerr << "Hmm read" << endl;

	out << hmm;
	if(out.good())
		cerr << "Hmm writtten" << endl;

	in.close();
	out.close();
	return 0;
}
