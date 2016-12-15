/*
 * PTU_test3.cpp
 *  test PTUnrooted evaluating at all posible nodes
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

using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 5) {
		cerr << "Evaluate a tree all possible nodes" << endl;
		cerr << "Usage:  " << argv[0] << " TREE-INFILE MSA-INFILE GTR-FILE OUTFILE" << endl;
		return -1;
	}

	ENABLE_INFO();

	ifstream in(argv[1]);
	ifstream msaIn(argv[2]);
	ifstream gtrIn(argv[3]);
	ofstream out(argv[4], ios_base::out | ios_base::binary);

	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	if(!msaIn.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}
	if(!gtrIn.is_open()) {
		cerr << "Unable to open " << argv[3] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[4] << endl;
		return -1;
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
		return -1;
	}
	else if(nRead != nLeaves) {
		cerr << "Loaded in " << nRead << " nodes from MSA but expecting " << nLeaves << " leaves in the PhyloTree " << endl;
		return -1;
	}
	else {
		infoLog << "MSA loaded successfully, read in " << nRead << " nodes with " << tree.numAlignSites() << " numSites" << endl;
	}

	GTR model;
	gtrIn >> model;
	infoLog << "DNA model loaded" << endl;
	infoLog << "MAX_COST_EXP: " << PhyloTreeUnrooted::MAX_COST_EXP << endl;

	clock_t t1 = clock();
	tree.initInCost();
	tree.initLeafCost(model);
	clock_t t2 = clock();
	infoLog << "Tree cost initiated, total eclipsed time: " << (t2 - t1) / static_cast<float>(CLOCKS_PER_SEC) << endl;

	for(size_t i = 0; i < tree.numNodes(); ++i) {
		tree.setRoot(i);
		infoLog << "Evaluating tree at node " << i << endl;
		tree.evaluate(model);
		double treeCost = tree.treeCost(model);
	}
	clock_t t3 = clock();
	infoLog << "All possible tree evaluated, total eclipsed time: " << (t3 - t1) / static_cast<float>(CLOCKS_PER_SEC) << endl;

	if(tree.save(out))
		infoLog << "Tree saved successfully" << endl;
	else {
		cerr << "Unable to save tree" << endl;
		return -1;
	}
}
