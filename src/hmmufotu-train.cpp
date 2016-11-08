/*
 * hmmufotu-train.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: zhengqi
 */


#include <iostream>
#include <fstream>
#include "HmmUFOtu.h"

using namespace std;
using namespace Eigen;
using namespace EGriceLab;
using namespace EGriceLab::Math;

static const double DEFAULT_SYMFRAC = 0.5;
static const string DEFAULT_DM_FILE = "gg_99_otus.dm";

/**
 * Print the usage information of this program
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <MSA-DB> [options]" << endl
		 << "Options:    -o FILE          : write output to FILE instead of stdout" << endl
		 << "            -symfrac DOUBLE  : conservation threshold for considering a site as a Match state in HMM [" << DEFAULT_SYMFRAC << "]" << endl
		 << "            -dm FILE         : use customized trained Dirichlet Model in FILE instead of build-in file " << endl
		 << "            -h|--help        : print this message and exit" << endl;
}

int main(int argc, char *argv[]) {
	ifstream in, dmin;
	ofstream of;
	double symfrac = DEFAULT_SYMFRAC;
	string infn;
	string outfn;
	string dmfn = DEFAULT_DM_FILE;
	string line;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printUsage(argv[0]);
		return -1;
	}
	if(cmdOpts.numMainOpts() != 1) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return -1;
	}

	infn = cmdOpts.getMainOpt(0);
	in.open(infn.c_str());
	if(!in.is_open()) {
		cerr << "Unable to open " << infn << endl;
		return -1;
	}

	if(cmdOpts.hasOpt("-o")) {
		outfn = cmdOpts.getOpt("-o");
		of.open(outfn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to " << outfn << endl;
			return -1;
		}
	}
	ostream& out = outfn.empty() ? cout : of;

	if(cmdOpts.hasOpt("-symfrac"))
		symfrac = atof(cmdOpts.getOpt("-symfrac").c_str());
	if(!(symfrac >= 0 && symfrac <= 1)) {
		cerr << "-symfrac must between 0 and 1" << endl;
	}

	if(cmdOpts.hasOpt("-dm"))
		dmfn = cmdOpts.getOpt("-dm");
	dmin.open(dmfn.c_str());
	if(!dmin.is_open()) {
		cerr << "Unable to open " << dmfn << endl;
		return -1;
	}

	/* Load in BandedHmmPrior for the HMM training */
	BandedHMMP7Prior hmmPrior;
	dmin >> hmmPrior;
	if(dmin.bad()) {
		cerr << "Failed to read in the HMM Prior file " << endl;
		cerr << hmmPrior << endl;
		return -1;
	}

	MSA* msa = MSA::load(in);
	if(!in.good()) {
		cerr << "Unable to load MSA database" << endl;
		return -1;
	}
	else
		cerr << "MSA loaded, found " << msa->getNumSeq() << " X " << msa->getCSLen() << " alignments" << endl;

	BandedHMMP7 hmm = BandedHMMP7::build(msa, symfrac, hmmPrior);

	cerr << "HMM trained" << endl;

	out << hmm;

	in.close();
	if(!outfn.empty())
	return 0;
}

