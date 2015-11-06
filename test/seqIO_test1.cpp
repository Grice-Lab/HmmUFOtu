#include <iostream>
#include <fstream>
#include "SeqIO.h"

using namespace std;
using namespace EGriceLab;
int main(int argc, char *argv[]) {
	if(argc != 5) {
		cerr << "Usage:  " << argv[0] << " INFILE OUTFILE INFORMAT OUTFORMAT" << endl;
		return -1;
	}

	SeqIO seqIn(argv[1], "dna", argv[3]);
	SeqIO seqOut(argv[2], "dna", argv[4], SeqIO::WRITE);

	while(seqIn.hasNext())
		seqOut.writeSeq(seqIn.nextSeq());

	return 0;
}
