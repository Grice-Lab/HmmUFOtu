/*
 * dna_model_test1.cpp
 *
 *  Created on: Sep 22, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include "MSA.h"
#include "PhyloTreeUnrooted.h"
#include "GTR.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 4) {
		cerr << "Usage:  " << argv[0] << " TREE-INFILE MSA-INFILE OUTFILE" << endl;
		return -1;
	}

	ifstream in(argv[1]);
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

	MSA msa;
	if(msa.load(msaIn))
		cerr << "MSA database loaded" << endl;
	else {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}

	EGriceLab::NT NTree;
	in >> NTree;
	cerr << "Newick Tree read" << endl;

	PTUnrooted tree(NTree);
	cerr << "PhyloTreeUnrooted constructed, total " << tree.numNodes() << " nodes found" << endl;

	long nLeaves = tree.numLeaves();
	long nRead = tree.loadMSA(msa);
	if(nRead == -1) {
		cerr << "Unable to read PhyloTree" << endl;
		return -1;
	}
	else if(nRead != nLeaves) {
		cerr << "Loaded in " << nRead << " nodes from MSA but expecting " << nLeaves << " leaves in the PhyloTree " << endl;
		return -1;
	}
	else {
		cerr << "MSA loaded into tree successfully, read in " << nRead << " nodes with " << tree.numAlignSites() << " numSites" << endl;
	}

	/* train an GTR model */
	GTR model;
	model.trainParams(tree.getModelTransitionSet(), tree.getModelFreqEst());

	cerr << "GTR model trained" << endl;
	out << model;
}
