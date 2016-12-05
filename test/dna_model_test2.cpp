/*
 * dna_model_test1.cpp
 *
 *  Created on: Sep 22, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include "MSA.h"
#include "PhyloTree.h"
#include "DNASubModel.h"
#include "GTR.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 5) {
		cerr << "Usage:  " << argv[0] << " TREE-INFILE MSA-INFILE DNA-MODELFILE OUTFILE" << endl;
		return -1;
	}

	ifstream msaIn(argv[2]);
	ifstream modelIn(argv[3]);
	ofstream out(argv[4]);
	if(!msaIn.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}
	if(!modelIn.is_open()) {
		cerr << "Unable to open " << argv[3] << endl;
		return -1;
	}

	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[3] << endl;
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

	/* read in DNA model */
	DNASubModel* model = new GTR();
	modelIn >> *model;

	if(modelIn.bad()) {
		cerr << "Failed to read in the DNA model" << endl;
		return -1;
	}
	else
		cerr << "DNA model read in successfully" << endl;

	cerr << "Starting evaluating tree" << endl;
	clock_t t1 = clock();

	model->evaluate(tree);

	clock_t t2 = clock();
	cerr << "Tree evaluated, time used: " << (t2 - t1) / CLOCKS_PER_SEC << endl;

	double c = model->cost(tree);
	cerr << "Tree cost calculated, c:" << c << endl;
}
