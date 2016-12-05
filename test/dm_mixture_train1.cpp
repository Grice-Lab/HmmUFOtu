#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <cstdlib>
#include "MSA.h"
#include "DirichletMixture.h"

using namespace std;
using namespace Eigen;
using namespace EGriceLab;
using namespace EGriceLab::Math;

int main(int argc, char *argv[]) {
	if(argc != 4) {
		cerr << "Usage:  " << argv[0] << " MSA-DB NUM-MIXTURE DM-OUTFILE" << endl;
		return -1;
	}

	ifstream in(argv[1]);
	int L = ::atoi(argv[2]);
	ofstream out(argv[3]);
	if(!in.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		return -1;
	}
	if(!out.is_open()) {
		cerr << "Unable to write to " << argv[2] << endl;
		return -1;
	}

	cerr << "File opened" << endl;
	MSA msa;
	if(msa.load(in))
		cerr << "MSA database loaded" << endl;
	else {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}

	MatrixXd data(msa.getAbc()->getSize(), msa.getCSLen());
	unsigned M = 0;
	for(unsigned j = 0; j < data.cols(); ++j) {
		if(msa.symWFrac(j) >= 0.5) // a match position conserved region
			data.col(M++) = msa.symWFreq(j);
	}
	data.conservativeResize(data.rows(), M);
	cerr << "Residual frequency data read, data:" << endl << data.rowwise().sum() << endl;

	DirichletMixture dm(data.rows(), L);
	cerr << "Dirichlet Mixture model constructed" << endl;
	double c = dm.trainML(data);
	cerr << "Dirichlet Mixture model EM trained" << endl;
	out << "Final cost: " << c << endl;
	out << dm << endl;
	cerr << "Dirichlet Mixture model written" << endl;

	return 0;
}
