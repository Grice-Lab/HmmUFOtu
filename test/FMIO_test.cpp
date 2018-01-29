#include <iostream>
#include <fstream>
#include "HmmUFOtuConst.h"
#include "HmmUFOtuEnv.h"
#include "ProgLog.h"
#include "CSFMIndex.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::HmmUFOtu;

int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " CSFM-INFILE CSFM-OUTFILE" << endl;
		return EXIT_FAILURE;
	}


	ifstream in(argv[1]);
	if(!in.is_open()) {
		cerr << "Unable to open CSFM file: " << argv[1] << endl;
		return EXIT_FAILURE;
	}

	ofstream out(argv[2]);
	if(!out.is_open()) {
		cerr << "Unable to write to CSFM file: " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	if(loadProgInfo(in).bad())
		return EXIT_FAILURE;
	CSFMIndex csfm;
	csfm.load(in);

	if(in.bad()) {
		cerr << "Unable to load CSFM file: " << argv[1] << endl;
	}
	else
		cerr << "CSFM index loaded" << endl;

	saveProgInfo(out);
	csfm.save(out);
	if(!out.good()) {
		cerr << "Unable to save CSFM index file: " << argv[2] << endl;
		return EXIT_FAILURE;
	}
	else
		cerr << "CSFM index file saved" << endl;

	return 0;
}
