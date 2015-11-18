#include <iostream>
#include <fstream>
#include "HmmUFOtuEnv.h"
#include "CSFMIndex.h"

using namespace std;
using namespace EGriceLab;
int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " CS-INFILE CS-OUTFILE" << endl;
		return -1;
	}

	ifstream in(argv[1]);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	ofstream out(argv[2]);
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}

	CSFMIndex* fmIdx = CSFMIndex::load(in);

	printVerboseInfo("FM index loaded");

	if(!fmIdx->save(out)) {
		cerr << "Unable to save CS-FMindex database" << endl;
		return -1;
	}

	return 0;
}
