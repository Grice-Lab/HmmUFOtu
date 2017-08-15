/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
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
static const string DEFAULT_SM_TYPE = "GTR";
static const string DEFAULT_TRAINING_METHOD = "Gojobori";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Train a DNA Substitution Model with customized data" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <MSA-FILE> <TREE-FILE> [options]" << endl
		 << "MSA-FILE  FILE                   : a multiple-alignment sequence file or pre-build MSA DB FILE" << endl
		 << "TREE-FILE  FILE                  : phylogenetic-tree file build on the MSA sequences" << endl
		 << "Options:    -o FILE              : write output to FILE instead of stdout" << endl
		 << "            --fmt  STR           : MSA format, supported format: 'fasta', 'msa'" << endl
		 << "            -s|--sub-model STR   : build a time-reversible DNA Substitution Model type, must be one of GTR, TN93, HKY85, F81, K80 or JC69 [" << DEFAULT_SM_TYPE << "]" << endl
		 << "            -m|--method  STR     : model training method using known phylogenetic tree data, either 'Gojobori' or 'Goldman' [" << DEFAULT_TRAINING_METHOD << "]" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string msaFn, treeFn, outFn;
	ifstream msaIn, treeIn;
	ofstream of;
	string fmt;
	string smType = DEFAULT_SM_TYPE;
	string method = DEFAULT_TRAINING_METHOD;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.empty() || cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printIntro();
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

	if(!(StringUtils::endsWith(treeFn, ".tree"))) {
		cerr << "Unrecognized TREE-FILE format, must be in Newick format" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-o"))
		outFn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("--fmt"))
		fmt = cmdOpts.getOpt("--fmt");

	if(cmdOpts.hasOpt("-s"))
		smType = cmdOpts.getOpt("-s");
	if(cmdOpts.hasOpt("--sub-model"))
		smType = cmdOpts.getOpt("--sub-model");

	if(cmdOpts.hasOpt("-m"))
		method = cmdOpts.getOpt("-m");
	if(cmdOpts.hasOpt("--method"))
		method = cmdOpts.getOpt("--method");

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* guess input format */
	if(fmt.empty()) {
		if(StringUtils::endsWith(msaFn, ".fasta") || StringUtils::endsWith(msaFn, ".fas")
		|| StringUtils::endsWith(msaFn, ".fa") || StringUtils::endsWith(msaFn, ".fna"))
			fmt = "fasta";
		else if(StringUtils::endsWith(msaFn, ".msa"))
			fmt = "msa";
		else {
			cerr << "Unrecognized format of MSA file '" << msaFn << "'" << endl;
			return EXIT_FAILURE;
		}
	}
	if(!(fmt == "fasta" || fmt == "msa")) {
		cerr << "Unsupported sequence format '" << fmt << "'" << endl;
		return EXIT_FAILURE;
	}

	/* open input files */
	if(fmt != "msa")
		msaIn.open(msaFn.c_str());
	else
		msaIn.open(msaFn.c_str(), ios_base::binary);
	if(!msaIn.is_open()) {
		cerr << "Unable to open '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

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
		msa.load(msaIn);
	}
	else if(msa.loadMSA(ALPHABET, msaIn, fmt) >= 0)
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
