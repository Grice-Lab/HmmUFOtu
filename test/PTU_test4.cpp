/*
 * PTU_test4.cpp
 *  PTUnrooted read placing test1
 *  Created on: Jan 5, 2017
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
	if(argc != 3) {
		cerr << "Placing the first leaf at all possible branches" << endl;
		cerr << "Usage:  " << argv[0] << " PTU-INFILE OUTFILE" << endl;
		return -1;
	}

	ENABLE_INFO();

	ifstream in(argv[1], ios_base::in | ios_base::binary);
	ofstream out(argv[2]);

	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to open " << argv[2] << endl;
		return -1;
	}

	PTUnrooted tree;
	tree.load(in);
	if(!in.bad()){
		cerr << "PTUnrooted loaded successfully, total " << tree.numNodes() << " nodes loaded" << endl;
		cerr << "Root id: " << tree.getRoot()->getId() << endl;
	}
	else {
		cerr << "Unable to load tree" << endl;
		return -1;
	}
	size_t N = tree.numNodes();

	infoLog << "Trying to place first leaf" << endl;
	DigitalSeq seq;
	for(int i = 0; i < N; ++i) {
		cerr << i << endl;
		const PTUnrooted::PTUNodePtr& node = tree.getNode(i);
		if(node->isLeaf()) {
			out << "Placing read copied from node " << i << " name: " << node->getName() << endl;
			seq = tree.getNode(i)->getSeq();
			break;
		}
	}

	for(int i = 0; i < N; ++i) {
		PTUnrooted::PTUNodePtr node = tree.getNode(i);
		if(node->isRoot())
			continue;
//		infoLog << "Copying subtree at node " << i << endl;
		PTUnrooted subtree = tree.copySubTree(node, node->getParent());
//		infoLog << "Placing read at subtree" << endl;
		double cost = subtree.placeSeq(seq, subtree.getNode(1), subtree.getNode(0), 0.03);
		out << "Read placed at node " << node->getId() << " name: " << node->getName()
				<< " cost: " << cost << endl;
	}
}
