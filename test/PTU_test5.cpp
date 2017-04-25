/*
 * PTU_test5.cpp
 *  PTUnrooted annotation test
 *  Created on: Jan 30, 2017
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <cfloat>
#include "MSA.h"
#include "PhyloTreeUnrooted.h"
#include "GTR.h"
#include "ProgLog.h"
#include "DNASubModelFactory.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 3) {
		cerr << "Annotate a PTUnrooted tree database" << endl;
		cerr << "Usage:  " << argv[0] << " PTU-INFILE PTU-OUTFILE" << endl;
		return EXIT_FAILURE;
	}

	INCREASE_LEVEL();

	ifstream in(argv[1], ios_base::in | ios_base::binary);
	ofstream out(argv[2], ios_base::out | ios_base::binary);

	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	PTUnrooted tree;
	tree.load(in);
	if(!in.bad())
		infoLog << "PTUnrooted loaded successfully, total " << tree.numNodes() << " nodes loaded" << endl;
	else {
		cerr << "Unable to load tree" << endl;
		return EXIT_FAILURE;
	}

	/* annotate tree */
	tree.annotate();
	infoLog << "PTUnrooted annotated" << endl;

	tree.save(out);
	if(!out.bad())
		infoLog << "PTUnrooted saved successfully" << endl;
	else {
		cerr << "Unable to save tree" << endl;
		return EXIT_FAILURE;
	}
}
