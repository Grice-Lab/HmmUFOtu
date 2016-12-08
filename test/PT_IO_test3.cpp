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
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " TREE-INFILE MSA-INFILE" << endl;
		return -1;
	}

	ifstream msaIn(argv[2]);
	if(!msaIn.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}

	MSA msa;
	if(msa.load(msaIn))
		cerr << "MSA database loaded" << endl;
	else {
		cerr << "Unable to load MSA database" << endl;
		return -1;
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
