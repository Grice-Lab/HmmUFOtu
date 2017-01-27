/*
 * PTU_test2.cpp
 *  test PTUnrooted evaluating at fixed root
 *  Created on: Dec 8, 2016
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
		cerr << "Evaluate a tree at given root" << endl;
		cerr << "Usage:  " << argv[0] << " TREE-INFILE MSA-INFILE OUTFILE" << endl;
		return -1;
	}

	ifstream in(argv[1]);
	ifstream msaIn(argv[2]);
	ofstream out(argv[3], ios_base::out | ios_base::binary);

	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	if(!msaIn.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[3] << endl;
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
	cerr << "NTRoot name: " << NTree.name << " neighbors: " << NTree.children.size() << endl;
	cerr << "PTRoot id: " << tree.getRoot()->getId() << " name: " << tree.getRoot()->getName()
			<< " neighbors: " << tree.getRoot()->numNeighbors() << endl;

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
		cerr << "MSA loaded successfully, read in " << nRead << " nodes with " << tree.numAlignSites() << " numSites" << endl;
	}

	ifstream modelIn("99_otus.gtr");
	GTR model;
	modelIn >> model;
	cerr << "DNA model loaded" << endl;
	tree.setModel(model);
	cerr << "DNA model set" << endl;

	clock_t t1 = clock();
	tree.initInLoglik();
	tree.initLeafLoglik();
	clock_t t2 = clock();
	cerr << "Tree loglik initiated, total eclipsed time: " << (t2 - t1) / static_cast<float>(CLOCKS_PER_SEC) << endl;

	tree.loglik(tree.getRoot());
	double treeLoglik = tree.treeLoglik();
	clock_t t3 = clock();
	cerr << "Tree evaluated, total eclipsed time: " << (t3 - t1) / static_cast<float>(CLOCKS_PER_SEC) << endl;
	cerr << "Final tree loglik: " << treeLoglik << endl;

	if(tree.save(out))
		cerr << "Tree saved successfully" << endl;
	else {
		cerr << "Unable to save tree" << endl;
		return -1;
	}
}
