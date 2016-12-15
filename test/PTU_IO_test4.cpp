/*
 * PTU_IO_test3.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include "MSA.h"
#include "PhyloTreeUnrooted.h"
#include "GTR.h"
using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 4) {
		cerr << "Load, repair and save a PTUnrooted binary file" << endl;
		cerr << "Usage:  " << argv[0] << " GTR-INFILE PTU-DB-IN PTU-DB-OUT" << endl;
		return -1;
	}

	ifstream modelIn(argv[1]);
	ifstream in(argv[2], ios_base::in | ios_base::binary);
	ofstream out(argv[3], ios_base::out | ios_base::binary);

	if(!modelIn.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[3] << endl;
		return -1;
	}

	GTR model;
	modelIn >> model;

	PTUnrooted tree;
	if(tree.load(in)) {
		cerr << "PTUnrooted loaded successfully, total " << tree.numNodes() << " nodes loaded" << endl;
		cerr << "Root id: " << tree.getRoot()->getId() << endl;
	}
	else {
		cerr << "Unable to load tree" << endl;
		return -1;
	}
	for(size_t i = 0; i < tree.numNodes(); ++i) {
		tree.setRoot(i);
		infoLog << "Evaluating tree at node " << i << endl;
		tree.evaluate(model);
		double treeCost = tree.treeCost(model);
	}
	cerr << "All root evaluated" << endl;

	if(tree.save(out))
		cerr << "PTUnrooted saved successfully" << endl;
	else {
		cerr << "Unable to save tree" << endl;
		return -1;
	}
}
