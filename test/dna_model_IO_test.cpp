/*
 * dna_model_IO_test.cpp
 *
 *  Created on: Jun 2, 2017
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_phylo.h"
#include "HmmUFOtuEnv.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::HmmUFOtu;

int main(int argc, char *argv[]) {
	if(argc < 3) {
		cerr << "Usage:  " << argv[0] << " DNA-SM-INFILE DNA-SM-OUTFILE [MODEL-TYPE]" << endl;
		return EXIT_FAILURE;
	}

	ifstream in(argv[1]);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	ofstream out(argv[2]);
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	/* guess DNA model type */
	char type[64];
	if(argc < 4) {
		cerr << "Guessing DNA substitution model type ... ";
		string line;
		while(std::getline(in, line)) {
			if(StringUtils::startsWith(line, "Type:")) {
				sscanf(line.c_str(), "Type: %s", type);
			}
		}
		cerr << type << endl;
		in.clear(); /* clear any error bit */
		in.seekg(0, ios_base::beg);
	}

	DNASubModel* model = DNASubModelFactory::createModel(type);
	in >> *model;
	if(in.bad()) {
		cerr << "Unable to read DNA model file: " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	else
		cerr << "DNA model read" << endl;

	out << *model;
	if(out.good())
		cerr << "DNA model written" << endl;
	else {
		cerr << "Unable to write to DNA model file: " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	return 0;
}
