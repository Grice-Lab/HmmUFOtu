/*
 * PTU_test3.cpp
 *  build a PTUnrooted tree index that pre-evaluating at all possible nodes
 *  Created on: Dec 8, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include "MSA.h"
#include "PhyloTreeUnrooted.h"
#include "GTR.h"
#include "ProgLog.h"
#include "DNASubModelFactory.h"

using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 6) {
		cerr << "Evaluate a tree all possible nodes" << endl;
		cerr << "Usage:  " << argv[0] << " TREE-INFILE ANNO-INFILE MSA-INFILE GTR-FILE OUTFILE" << endl;
		return -1;
	}

	INCREASE_LEVEL();

	ifstream in(argv[1]);
	ifstream annoIn(argv[2]);
	ifstream msaIn(argv[3]);
	ifstream gtrIn(argv[4]);
	ofstream out(argv[5], ios_base::out | ios_base::binary);

	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	if(!annoIn.is_open()) {
		cerr << "Unable to open " << argv[2]<< endl;
		return EXIT_FAILURE;
	}
	if(!msaIn.is_open()) {
		cerr << "Unable to open " << argv[3] << endl;
		return EXIT_FAILURE;
	}
	if(!gtrIn.is_open()) {
		cerr << "Unable to open " << argv[4] << endl;
		return EXIT_FAILURE;
	}
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[5] << endl;
		return EXIT_FAILURE;
	}

	MSA msa;
	if(msa.load(msaIn))
		infoLog << "MSA database loaded" << endl;
	else {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}

	EGriceLab::NT NTree;
	in >> NTree;
	infoLog << "Newick Tree read" << endl;

	PTUnrooted tree(NTree);
	infoLog << "PhyloTreeUnrooted constructed, total " << tree.numNodes() << " nodes found" << endl;
	infoLog << "NTRoot name: " << NTree.name << " neighbors: " << NTree.children.size() << endl;
	infoLog << "PTRoot id: " << tree.getRoot()->getId() << " name: " << tree.getRoot()->getName()
			<< " neighbors: " << tree.getRoot()->numNeighbors() << endl;

	long nLeaves = tree.numLeaves();
	long nRead = tree.loadMSA(msa);
	if(nRead == -1) {
		cerr << "Unable to read PhyloTree" << endl;
		return EXIT_FAILURE;
	}
	else if(nRead != nLeaves) {
		cerr << "Loaded in " << nRead << " nodes from MSA but expecting " << nLeaves << " leaves in the PhyloTree " << endl;
		return EXIT_FAILURE;
	}
	else {
		infoLog << "MSA loaded successfully, read in " << nRead << " nodes with " << tree.numAlignSites() << " numSites" << endl;
	}

	tree.loadAnnotation(annoIn);
	if(annoIn.bad()) {
		cerr << "Failed to load annotation " << argv[2] << endl;
		return EXIT_FAILURE;
	}
	else
		cerr << "Taxa annotation loaded" << endl;

	tree.formatName();
	cerr << "Tree node names formatted" << endl;

	DNASubModel* model = DNASubModelFactory::createModel("GTR");
	gtrIn >> *model;
	infoLog << "DNA model loaded" << endl;
	tree.setModel(model);
	infoLog << "DNA model set" << endl;

	clock_t t1 = clock();
	tree.initBranchLoglik();
	tree.initLeafMat();
	clock_t t2 = clock();
	infoLog << "Tree loglik initiated, total eclipsed time: " << (t2 - t1) / static_cast<float>(CLOCKS_PER_SEC) << endl;

	EGriceLab::PTUnrooted::PTUNodePtr oldRoot = tree.getRoot();
	for(size_t i = 0; i < tree.numNodes(); ++i) {
		tree.setRoot(i);
		infoLog << "Evaluating tree at node " << i << endl;
		tree.evaluate();
	}
	/* reset root to original */
	tree.setRoot(oldRoot);
	clock_t t3 = clock();
	infoLog << "All possible tree evaluated, total eclipsed time: " << (t3 - t1) / static_cast<float>(CLOCKS_PER_SEC) << endl;

	if(tree.save(out))
		infoLog << "Tree saved successfully" << endl;
	else {
		cerr << "Unable to save tree" << endl;
		return 1;
	}
}
