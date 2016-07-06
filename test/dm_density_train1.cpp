#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "MSA.h"
#include "DirichletDensity.h"

using namespace std;
using namespace Eigen;
using namespace EGriceLab;
using namespace EGriceLab::Math;

int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "Usage:  " << argv[0] << " MSA-DB DM-OUTFILE" << endl;
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

	cerr << "File opened" << endl;
	MSA* msa = MSA::load(in);
	if(!in.good()) {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}
	else {
		cerr << "MSA database loaded" << endl;
	}

	MatrixXd data(msa->getAbc()->getSize(), msa->getCSLen());
	unsigned M = 0;
	for(unsigned j = 0; j < data.cols(); ++j) {
		if(msa->symWFrac(j) < 0.5) // an insertion, non-conserved region
			data.col(M++) = msa->symWFreq(j);
	}
	data.conservativeResize(data.rows(), M);
	cerr << "Residual frequency data read, total M:" << M << endl;

	DirichletDensity dm(data.rows(), 100);
	cerr << "Dirichlet Density model constructed" << endl;
	dm.trainML(data);
	cerr << "Dirichlet Density model Maximum Likelihood trained" << endl;
	out << dm << endl;
	cerr << "Dirichlet Density model written" << endl;

	return 0;
}
