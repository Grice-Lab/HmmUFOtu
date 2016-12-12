#include <iostream>
#include <fstream>
#include "MSA.h"

using namespace std;
using namespace EGriceLab;
int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " DB-INFILE DB-OUTFILE" << endl;
		return -1;
	}

	ifstream in(argv[1], ios_base::in | ios_base::binary);
	ofstream out(argv[2], ios_base::out | ios_base::binary);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return -1;
	}
	cerr << "File opened" << endl;

	MSA msa;
	if(msa.load(in))
		cerr << "MSA database loaded" << endl;
	else {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}

	cerr << "Total seqNum:" << msa.getNumSeq() << endl;

	if(!msa.save(out)) {
		cerr << "Unable to save MSA database" << endl;
		return -1;
	}
	else {
		cerr << "MSA database saved" << endl;
	}

	return 0;
}
