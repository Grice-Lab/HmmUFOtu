#include <iostream>
#include <fstream>
#include "MSA.h"

using namespace std;
using namespace EGriceLab;
int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " DB-INFILE OUTFILE" << endl;
		return -1;
	}

	ifstream in(argv[1]);
	ofstream out(argv[2]);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return -1;
	}

	MSA* msa = MSA::load(in);
	if(!in.good()) {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}
	msa->prune();
	cerr << "MSA pruned" << endl;

	cerr << "Effective seqNum:" << msa->getEffectSeqNum() << endl;

	out << "pos\tidentity\tsym_frac\tgap_frac" << endl;
	for(unsigned j = 0; j != msa->getCSLen(); ++j)
		out << j << "\t" << msa->identityAt(j) << "\t" <<
			msa->symFrac(j) << "\t" << msa->gapFrac(j) << endl;

	return 0;
}
