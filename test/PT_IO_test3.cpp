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

	clock_t t1 = clock();
	EGriceLab::PT tree;
	int nRead = tree.readTree(argv[1], "newick", msa);
	if(nRead != -1)
		cerr << "PhyloTree read with " << nRead << " assigned seq" << endl;
	else {
		cerr << "Unable to read PhyloTree" << endl;
		return -1;
	}
	clock_t t2 = clock();
	cerr << "Time elapsed " << (t2 - t1) / static_cast<float>(CLOCKS_PER_SEC) << endl;

//	tree.setIDandParent();
//	tree.annotate();

	clock_t t3 = clock();
	cerr << "PhyloTree annotated" << endl;
	cerr << "Time elapsed " << (t3 - t2) / static_cast<float>(CLOCKS_PER_SEC) << endl;

	PhyloTree tree2(tree);
	tree2.updateParent();
	clock_t t4 = clock();
	cerr << "PhyloTree copied" << endl;
	cerr << "Time elapsed " << (t4 - t3) / static_cast<float>(CLOCKS_PER_SEC) << endl;

	out << tree << endl;
}
