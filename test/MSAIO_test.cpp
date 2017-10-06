#include <iostream>
#include <fstream>
#include "MSA.h"
#include "HmmUFOtuConst.h"
#include "HmmUFOtuEnv.h"

using namespace std;
using namespace EGriceLab;
int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " DB-INFILE DB-OUTFILE" << endl;
		return EXIT_FAILURE;
	}

	INCREASE_LEVEL();

	ifstream in(argv[1], ios_base::in | ios_base::binary);
	ofstream out(argv[2], ios_base::out | ios_base::binary);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return EXIT_FAILURE;
	}
	cerr << "File opened" << endl;

	if(loadProgInfo(in).bad())
		return EXIT_FAILURE;
	MSA msa;
	msa.load(in);
	if(!in.bad())
		cerr << "MSA database loaded" << endl;
	else {
		cerr << "Unable to load MSA database: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	cerr << "Total seqNum:" << msa.getNumSeq() << endl;

	saveProgInfo(out);
	msa.save(out);
	if(out.bad()) {
		cerr << "Unable to save MSA database: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	else
		cerr << "MSA database saved" << endl;

	return 0;
}
