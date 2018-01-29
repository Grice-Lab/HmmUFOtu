/*
 * PTU_IO_test.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include "MSA.h"
#include "PhyloTreeUnrooted.h"
#include "HmmUFOtuConst.h"
#include "HmmUFOtuEnv.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::HmmUFOtu;

int main(int argc, const char* argv[]) {
	if(argc != 3) {
		cerr << "Load and save a PTUnrooted binary file" << endl;
		cerr << "Usage:  " << argv[0] << " PTU-DB-IN PTU-DB-OUT" << endl;
		return EXIT_FAILURE;
	}

	ifstream in(argv[1], ios_base::in | ios_base::binary);
	ofstream out(argv[2], ios_base::out | ios_base::binary);

	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	if(loadProgInfo(in).bad())
		return EXIT_FAILURE;
	PTUnrooted tree;
	tree.load(in);
	if(!in.bad()){
		cerr << "PTUnrooted loaded successfully, total " << tree.numNodes() << " nodes loaded" << endl;
		cerr << "Root id: " << tree.getRoot()->getId() << endl;
	}
	else {
		cerr << "Unable to load tree" << endl;
		return EXIT_FAILURE;
	}

	saveProgInfo(out);
	tree.save(out);

	if(!out.bad())
		cerr << "PTUnrooted saved successfully" << endl;
	else {
		cerr << "Unable to save tree" << endl;
		return EXIT_FAILURE;
	}
}
