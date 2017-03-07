/*
 * hmmufotu-train-sm.cpp
 *  train a customized DNA Substitution Model using a phylogenetic tree and associated MSA database
 *
 *  Created on: Feb 9, 2017
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <string>
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_phylo.h"
#include "EGMath.h"

using namespace std;
using namespace EGriceLab;

/** default values */
static const string ALPHABET = "dna";
static const string DEFAULT_SUB_MODEL_TYPE = "GTR";
static const string DEFAULT_TRAINING_METHOD = "Gojobori";

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Train a customized DNA Substitution Model for " << progName << " analysis" << endl
		 << "Usage:    " << progName << "  <MSA-FILE> <TREE-FILE> [options]" << endl
		 << "MSA-FILE  FILE                   : a multiple-alignment sequence file or pre-build MSA DB FILE" << endl
		 << "TREE-FILE  FILE                  : phylogenetic-tree file build on the MSA sequences" << endl
		 << "Options:    -o FILE              : write output to FILE instead of stdout" << endl
		 << "            -s|--sub-model STR   : build a time-reversible DNA Substitution Model type, must be one of GTR, TN93, HKY85, F81, K80 or JC69 [" << DEFAULT_SUB_MODEL_TYPE << "]" << endl
		 << "            -m|--method  STR     : model training method using known phylogenetic tree data, either 'Gojobori' or 'Goldman' [" << DEFAULT_TRAINING_METHOD << "]" << endl
		 << "            -v  FLAG             : enable verbose information" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string msaFn, treeFn, outFn;
	ifstream msaIn, treeIn;
	ofstream of;
	string fmt;
	string smType = DEFAULT_SUB_MODEL_TYPE;
	string method = DEFAULT_TRAINING_METHOD;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(cmdOpts.numMainOpts() != 2) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	msaFn = cmdOpts.getMainOpt(0);
	treeFn = cmdOpts.getMainOpt(1);

	/* check input format */
	if(StringUtils::endsWith(msaFn, ".fasta") || StringUtils::endsWith(msaFn, ".fas")
		|| StringUtils::endsWith(msaFn, ".fa") || StringUtils::endsWith(msaFn, ".fna"))
		fmt = "fasta";
	else if(StringUtils::endsWith(msaFn, ".msa"))
		fmt = "msa";
	else {
		cerr << "Unrecognized format of MSA file '" << msaFn << "'" << endl;
		return EXIT_FAILURE;
	}

	if(!(StringUtils::endsWith(treeFn, ".tree"))) {
		cerr << "Unrecognized TREE-FILE format, must be in Newick format" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-o"))
		outFn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("-s"))
		smType = cmdOpts.getOpt("-s");
	if(cmdOpts.hasOpt("--sub-model"))
		smType = cmdOpts.getOpt("--sub-model");

	if(cmdOpts.hasOpt("-m"))
		method = cmdOpts.getOpt("-m");
	if(cmdOpts.hasOpt("--method"))
		method = cmdOpts.getOpt("--method");

	if(cmdOpts.hasOpt("-v"))
		ENABLE_INFO();

	/* open input files */
	treeIn.open(treeFn.c_str());
	if(!treeIn.is_open()) {
		cerr << "Unable to open '" << treeFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* open output files */
	if(!outFn.empty()) {
		of.open(outFn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to '" << outFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}
	ostream& out = of.is_open() ? of : cout;

	/* Load data */
	MSA msa;
	if(fmt == "msa") { /* binary file provided */
		ifstream in(msaFn.c_str());
		msa.load(in);
		if(!in.good()) {
			cerr << "Unable to load MSA database from '" << msaFn << "'" << endl;
			return EXIT_FAILURE;
		}
	}
	else if(msa.loadMSAFile(ALPHABET, msaFn, fmt) >= 0)
		infoLog << "MSA loaded" << endl;
	else {
		cerr << "Unable to load MSA seq from '" << msaFn << "'" << endl;
		return EXIT_FAILURE;
	}

	if(!msa.pruned()) {
		msa.prune(); /* prune MSA if necessary*/
		infoLog << "MSA pruned" << endl;
	}
	infoLog << "MSA database created for " << msa.getNumSeq() << " X " << msa.getCSLen() << " aligned sequences" << endl;

	/* build ptu */
	EGriceLab::NewickTree NTree;
	treeIn >> NTree;
	if(treeIn.bad()) {
		cerr << "Unable to read Newick tree in '" << treeFn << "'" << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Newick Tree read" << endl;

	PTUnrooted tree(NTree);
	infoLog << "Phylogenetic Tree constructed with total " << tree.numNodes() << " nodes" << endl;

	size_t nLeaves = tree.numLeaves();
	size_t nRead = tree.loadMSA(msa);
	if(nRead == -1) {
		cerr << "Unable to load MSA into Phylogenetic Tree" << endl;
		return EXIT_FAILURE;
	}
	else if(nRead != nLeaves) {
		cerr << "Unmatched MSA and Tree. Found " << nRead << " sequences from MSA but expecting " << nLeaves << " leaves in the Phylogenetic Tree " << endl;
		return EXIT_FAILURE;
	}
	else
		infoLog << "MSA loaded into Phylogenetic Tree" << endl;

	/* train DNA sub model */
	DNASubModel* model = DNASubModelFactory::createModel(smType);
	model->trainParams(tree.getModelTransitionSet(method), tree.getModelFreqEst());
	infoLog << "DNA Substitution Model trained" << endl;

	/* output */
	if(out << *model)
		infoLog << "Model written" << endl;
	else {
		cerr << "Unable to write model: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
}
