/*
 * PTU_IO_test4.cpp
 *  Load PTUnrooted and test the tree costs
 *  Created on: Dec 16, 2016
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include "MSA.h"
#include "PhyloTreeUnrooted.h"
#include "GTR.h"
using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 3) {
		cerr << "Test tree cost on every node on a saved PTUnrooted tree" << endl;
		cerr << "Usage:  " << argv[0] << " PTU-DB-IN COST-OUT" << endl;
		return 1;
	}

	ifstream in(argv[1], ios_base::in | ios_base::binary);
	ofstream out(argv[2]);

	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return EXIT_FAILURE;
	}

	PTUnrooted tree;
	tree.load(in);
	if(!in.bad()) {
		cerr << "PTUnrooted loaded successfully, total " << tree.numNodes() << " nodes loaded" << endl;
		cerr << "Root id: " << tree.getRoot()->getId() << endl;
	}
	else {
		cerr << "Unable to load tree" << endl;
		return EXIT_FAILURE;
	}

	double maxCostDiff = 0;
	double prevCost = 0;
	for(size_t i = 0; i < tree.numNodes(); ++i) {
		tree.setRoot(i);
		double cost = tree.treeCost();
		out << "Tree cost at node " << i << " isLeaf? " << tree.getRoot()->isLeaf() << " cost: " << cost << endl;
		if(prevCost != 0 && ::abs(cost - prevCost) > maxCostDiff)
			maxCostDiff = ::abs(cost - prevCost);
		prevCost = cost;
	}
	out << "Max cost difference of all nodes: " << maxCostDiff << endl;
	out << "Max relative cost difference of all nodes: " << maxCostDiff / prevCost << endl;
}
