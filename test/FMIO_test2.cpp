#include <iostream>
#include <fstream>
#include "HmmUFOtuEnv.h"
#include "ProgLog.h"
#include "CSFMIndex.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " CS-INFILE CS-OUTFILE" << endl;
		return -1;
	}
	VERBOSE_LEVEL += 2;

	ifstream in(argv[1]);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	else
		infoLog << "infile opened successfully" << endl;

	ofstream out(argv[2]);
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}
	else
		infoLog << "outfile opened successfully" << endl;

	CSFMIndex* fmIdx = CSFMIndex::load(in);

	infoLog << "FM index loaded" << endl;

	if(!fmIdx->save(out)) {
		cerr << "Unable to save CS-FMindex database" << endl;
		return -1;
	}
	else
		infoLog << "FM index saved successfully" << endl;

	return 0;
}
