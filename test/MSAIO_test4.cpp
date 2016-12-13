#include <iostream>
#include <fstream>
#include "MSA.h"
#include "SeqIO.h"

using namespace std;
using namespace EGriceLab;
int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " DB-INFILE FASTA-OUTFILE" << endl;
		return -1;
	}

	ifstream in(argv[1], ios_base::in | ios_base::binary);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
//	SeqIO out(argv[2], "dna", "fasta", SeqIO::WRITE);
	ofstream out(argv[2]);
	cerr << "File opened" << endl;

	MSA msa;
	if(msa.load(in))
		cerr << "MSA database loaded" << endl;
	else {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}

	cerr << "csLen: " << msa.getCSLen() << endl;
	cerr << "Total seqNum: " << msa.getNumSeq() << endl;

	for(unsigned i = 0; i < msa.getNumSeq(); ++i)
//		out.writeSeq(msa.primarySeqAt(i));
//		out << '>' << msa.seqNameAt(i) << endl << msa.primarySeqAt(i).removeGaps().getSeq() << endl;
		out << '>' << msa.seqNameAt(i) << endl << msa.seqAt(i) << endl;


	return 0;
}
