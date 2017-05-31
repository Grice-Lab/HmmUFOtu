#include <iostream>
#include <fstream>
#include "DirichletDensity.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::Math;
int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " DM-INFILE DM-OUTFILE" << endl;
		return EXIT_FAILURE;
	}

	ifstream in(argv[1]);
	ofstream out(argv[2]);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	DirichletDensity dd;
	in >> dd;

	if(!in.good()) {
		cerr << "Unable to read DM density file: " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	else
		cerr << "DM density file read" << endl;

	out << dd;
	if(!out.good()) {
		cerr << "Unable to write DM density file: " << argv[2] << endl;
		return EXIT_FAILURE;
	}
	else
		cerr << "DM density written" << endl;

	return 0;
}
