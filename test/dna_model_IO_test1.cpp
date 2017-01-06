/*
 * dna_model_IO_test1.cpp
 *
 *  Created on: Sep 22, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include "MSA.h"
#include "DNASubModel.h"
#include "GTR.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " INFILE OUTFILE" << endl;
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

	DNASubModel* model = new GTR();

	in >> *model;
	if(in.bad()) {
		cerr << "Failed to read in DNA model" << endl;
		return -1;
	}

	out << *model;
	if(out.fail()) {
		cerr << "Failed to write DNA model" << endl;
		return -1;
	}

}
