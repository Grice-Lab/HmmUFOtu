#include <iostream>
#include <fstream>
#include <string>
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

	cerr << "File opened" << endl;
	MSA* msa = MSA::load(in);
	if(!in.good()) {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}
	else {
		cerr << "MSA database loaded" << endl;
	}

	cerr << "Total numSeq:" << msa->getNumSeq() << endl;

	out << "pos";
	string sym = msa->getAbc()->getSymbol();
	for(string::const_iterator it = sym.begin(); it != sym.end(); ++it)
		out << " " << *it;
	out << "\n";
	unsigned L = msa->getCSLen();
	for(unsigned j = 0; j < L; ++j) {
		out << j << " " << msa->symFreq(j).transpose() << endl;
	}

	return 0;
}
