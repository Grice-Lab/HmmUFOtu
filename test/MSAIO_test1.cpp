#include <iostream>
#include <fstream>
#include "MSA.h"

using namespace std;
using namespace EGriceLab;
int main(int argc, char *argv[]) {
	if(argc != 4) {
		cerr << "Usage:  " << argv[0] << " MSA-INFILE format DB-OUTFILE" << endl;
		return -1;
	}

	ofstream out(argv[3], ios_base::out | ios_base::binary);
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[3] << endl;
		return -1;
	}

	MSA msa;
	if(msa.loadMSAFile("dna", argv[1], argv[2]) >= 0)
		cerr << "MSA file loaded" << endl;
	else {
		cerr << "Unable to load MSA file" << endl;
		return -1;
	}

	msa.prune();
	cerr << "MSA pruned" << endl;

	cerr << "csLen: " << msa.getCSLen() << endl;
	cerr << "Total seqNum: " << msa.getNumSeq() << endl;

	if(!msa.save(out)) {
		cerr << "Unable to save MSA database" << endl;
		return -1;
	}
	cerr << "MSA saved" << endl;

	return 0;
}
