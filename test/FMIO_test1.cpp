#include <iostream>
#include <fstream>
#include "CSFMIndex.h"

using namespace std;
using namespace EGriceLab;
int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " MSA-INFILE CS-OUTFILE" << endl;
		return -1;
	}

	string inFn(argv[1]);
	string outFn(argv[2]);

	ofstream out(outFn.c_str(), ios_base::out | ios_base::binary);
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[3] << endl;
		return EXIT_FAILURE;
	}

	string fmt;
	/* guess input format */
	if(StringUtils::endsWith(inFn, ".fasta") || StringUtils::endsWith(inFn, ".fas")
		|| StringUtils::endsWith(inFn, ".fa") || StringUtils::endsWith(inFn, ".fna"))
		fmt = "fasta";
	else if(StringUtils::endsWith(inFn, ".msa"))
		fmt = "msa";
	else {
		cerr << "Unrecognized MSA file format" << endl;
		return EXIT_FAILURE;
	}

	/* Load data */
	MSA msa;
	if(fmt == "msa") { /* binary file provided */
		ifstream in(inFn.c_str());
		msa.load(in);
		if(!in.good()) {
			cerr << "Unable to load MSA database from '" << inFn << "'" << endl;
			return EXIT_FAILURE;
		}
	}
	else if(msa.loadMSAFile("dna", inFn, fmt) < 0) {
		cerr << "Unable to load MSA seq from '" << inFn << "'" << endl;
		return EXIT_FAILURE;
	}

	msa.prune(); /* prune MSA if necessary*/
//	infoLog << "MSA database created for " << msa.getNumSeq() << " X " << msa.getCSLen() << " aligned sequences" << endl;

	CSFMIndex csfm;
	csfm.build(msa);

	cerr << "FM-index built" << endl;

	if(!csfm.save(out)) {
		cerr << "Unable to save CS-FMindex database" << endl;
		return EXIT_FAILURE;
	}

	return 0;
}
