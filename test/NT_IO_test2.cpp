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
#include "NewickTree.h"
using namespace std;
using namespace EGriceLab;

int main(int argc, const char* argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " TREE-INFILE OUTFILE" << endl;
		return -1;
	}

	ifstream in(argv[1]);
	ofstream out(argv[2]);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return -1;
	}

	EGriceLab::NT tree;
	in >> tree;

	out << tree;
}
