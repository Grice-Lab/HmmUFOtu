/*
 * hmmufotu-inspect.cpp
 * Inspect a hmmufotu database
 *  Created on: Feb 28, 2017
 *      Author: zhengqi
 */

#include <iostream>
#include "HmmUFOtu.h"

using namespace std;
using namespace EGriceLab;

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Inspect a HmmUFOtu database" << endl
		 << "Usage:    " << progName << "  <DB-NAME> [options]" << endl
		 << "DB-NAME  STR                    : HmmUFOtu database name (prefix)" << endl
		 << "Options:    -sm  FLAG           : report the embedded build-in or customized DNA Submission Model in database" << endl
		 << "            -dg  FLAG           : report the embedded build-in or customized Discrete Gamma Model (if enabled during training) in database" << endl
		 << "            -v  FLAG            : enable verbose information" << endl
		 << "            -h|--help           : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, msaFn, csfmFn, hmmFn, ptuFn;
	ifstream msaIn, csfmIn, hmmIn, ptuIn;
	bool showSm, showDg;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(cmdOpts.numMainOpts() != 1) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	dbName = cmdOpts.getMainOpt(0);
	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	if(cmdOpts.hasOpt("-sm"))
		showSm = true;

	if(cmdOpts.hasOpt("-dg"))
		showDg = true;

	msaFn = dbName + MSA_FILE_SUFFIX;
	csfmFn = dbName + CSFM_FILE_SUFFIX;
	hmmFn = dbName + HMM_FILE_SUFFIX;
	ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;

	/* open inputs */
	msaIn.open(msaFn.c_str(), ios_base::in | ios_base::binary);
	if(!msaIn) {
		cerr << "Unable to open MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	csfmIn.open(csfmFn.c_str(), ios_base::in | ios_base::binary);
	if(!csfmIn) {
		cerr << "Unable to open CSFM-index '" << csfmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	hmmIn.open(hmmFn.c_str());
	if(!hmmIn) {
		cerr << "Unable to open HMM profile '" << hmmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	ptuIn.open(ptuFn.c_str(), ios_base::in | ios_base::binary);
	if(!ptuIn) {
		cerr << "Unable to open PTU data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	int csLen;

	infoLog << "Inspecting MSA data ..." << endl;
	MSA msa;
	msa.load(msaIn);
	if(msaIn.bad()) {
		cerr << "Failed to load MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	csLen = msa.getCSLen();
	cout << "MSA loaded. Number of seq: " << msa.getNumSeq() << " CS length: " << csLen << endl;

	infoLog << "Inspecting CSFM-index ..." << endl;
	CSFMIndex csfm;
	csfm.load(csfmIn);
	if(csfmIn.bad()) {
		cerr << "Failed to load CSFM-index '" << csfmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	cout << "CSFM-index loaded. Concatednated length: " << csfm.getConcatLen() << " CS length: " << csfm.getCSLen() << endl;
	if(csfm.getCSLen() != csLen) {
		cerr << "Error: Unmatched CS length between CSFM-index and MSA data" << endl;
		return EXIT_FAILURE;
	}

	infoLog << "Inspecting HMM profile ..." << endl;
	BandedHMMP7 hmm;
	hmmIn >> hmm;
	if(hmmIn.bad()) {
		cerr << "Unable to read HMM profile '" << hmmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	cout << "HMM profile read. Name: " << hmm.getName() << " Alphabet: "
		 << hmm.getNuclAbc()->getAlias() << " Profile size: " << hmm.getProfileSize() << endl;
	if(hmm.getProfileSize() > csLen) {
		cerr << "Error: HMM profile size is found greater than the MSA CS length" << endl;
		return EXIT_FAILURE;
	}

	infoLog << "Inspecting Phylogenetic tree data ..." << endl;
	PTUnrooted ptu;
	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Unable to load Phylogenetic tree data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	cout << "Phylogenetic tree loaded. Root ID: " << ptu.getRoot()->getId()
		 << " Number of leaves: " << ptu.numLeaves()
		 << " Number of nodes: " << ptu.numNodes()
		 << " Number of branches: " << ptu.numBranches()
		 << " Number of sites: " << ptu.numAlignSites() << endl;
	cout << "Root name: " << ptu.getRoot()->getAnno() << endl
		 << "Entire tree log-liklihood: " << ptu.treeLoglik() << endl;

	if(showSm)
		cout << (*ptu.getModel());

	if(showDg && ptu.getDGModel() != NULL)
		cout << "Discrete Gamma Model is enabled for this tree" << endl
		     << "Number of categories used: " << ptu.getDGModel()->getK()
			 << "Shape parameter: " << ptu.getDGModel()->getShape() << endl;

}

