/*
 * bHMM_train_test1.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: zhengqi
 */


#include <iostream>
#include <fstream>
#include "HmmUFOtu.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " MSA-DBFILE HMM-OUTFILE" << endl;
		return 0;
	}

	ifstream in(argv[1]);
	ofstream out(argv[2]);

	MSA* msa = MSA::load(in);
	if(!in.good()) {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}
	cerr << msa->getNumSeq() << " x " << msa->getCSLen() << endl;

	BandedHMMP7 hmm = BandedHMMP7::build(msa, 0.5);

	cerr << "HMM trained" << endl;

	out << hmm;

	in.close();
	out.close();
	return 0;
}

