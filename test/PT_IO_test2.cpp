/*
 * PT_IO_test1.cpp
 *
 *  Created on: Aug 25, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include "MSA.h"
#include "PhyloTree.h"
using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 4) {
		cerr << "Usage:  " << argv[0] << " TREE-INFILE MSA-INFILE OUTFILE" << endl;
		return -1;
	}

	ifstream msaIn(argv[2]);
	ofstream out(argv[3]);
	if(!msaIn.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[3] << endl;
		return -1;
	}

	MSA* msa = MSA::load(msaIn);
	if(!msaIn.good()) {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}
	else {
		cerr << "MSA database loaded" << endl;
	}

	EGriceLab::PT tree;
	int nRead = tree.readTree(argv[1], "newick", msa);
	if(nRead != -1)
		cerr << "Read in PhyloTree with " << nRead << " assigned seq" << endl;
	else {
		cerr << "Unable to read PhyloTree" << endl;
		return -1;
	}

	out << tree << endl;
}
