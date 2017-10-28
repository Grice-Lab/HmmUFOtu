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
 * hmmufotu-inspect.cpp
 * Inspect a hmmufotu database
 *  Created on: Feb 28, 2017
 *      Author: zhengqi
 */

#include <iostream>
#include <boost/lexical_cast.hpp>
#include "HmmUFOtu.h"

using namespace std;
using namespace EGriceLab;

static const string TREE_FORMAT = "newick";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Inspect an HmmUFOtu database, and optionally export its contents" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <DB-NAME> [options]" << endl
		 << "DB-NAME  STR                    : HmmUFOtu database name (prefix)" << endl
		 << "Options:    -sm  FLAG           : report the embedded build-in or customized DNA Submission Model in database" << endl
		 << "            -dg  FLAG           : report the embedded build-in or customized Discrete Gamma Model (if enabled during training) in database" << endl
		 << "            -t|--tree  FILE     : write the phylogenetic tree of this database to FILE in Newick format" << endl
		 << "            -a|--anno  FILE     : write the tree node taxonomy annoation of this database to FILE" << endl
		 << "            -s|--seq  FILE      : write the multiple-sequence alignment of this database to FILE in fasta format" << endl
		 << "            --use-dbname  FLAG  : use DBNAME as prefix for all tree nodes" << endl
		 << "            -n|--node  FLAG     : write sequence alignment of all nodes instead of just leaves, ignored if -s is not set" << endl
		 << "            -v  FLAG            : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version          : show program version and exit" << endl
		 << "            -h|--help           : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, msaFn, csfmFn, hmmFn, ptuFn;
	string treeFn, annoFn, seqFn;
	ifstream msaIn, csfmIn, hmmIn, ptuIn;
	ofstream treeOut, annoOut, seqOut;
	SeqIO seqO;
	bool showSm = false;
	bool showDg = false;
	bool leafOnly = true;
	bool useDBName = false;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.empty() || cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printIntro();
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(cmdOpts.hasOpt("--version")) {
		printVersion(argv[0]);
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

	if(cmdOpts.hasOpt("-t"))
		treeFn = cmdOpts.getOpt("-t");
	if(cmdOpts.hasOpt("--tree"))
		treeFn = cmdOpts.getOpt("--tree");

	if(cmdOpts.hasOpt("-a"))
		annoFn = cmdOpts.getOpt("-a");
	if(cmdOpts.hasOpt("--anno"))
		annoFn = cmdOpts.getOpt("--anno");

	if(cmdOpts.hasOpt("-s"))
		seqFn = cmdOpts.getOpt("-s");
	if(cmdOpts.hasOpt("--seq"))
		seqFn = cmdOpts.getOpt("--seq");

	if(cmdOpts.hasOpt("--use-dbname"))
		useDBName = true;

	if(cmdOpts.hasOpt("-n") || cmdOpts.hasOpt("--node"))
		leafOnly = false;

	msaFn = dbName + MSA_FILE_SUFFIX;
	csfmFn = dbName + CSFM_FILE_SUFFIX;
	hmmFn = dbName + HMM_FILE_SUFFIX;
	ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;
	string nodePrefix = !useDBName ? "" : dbName + "_";

	/* open inputs */
	msaIn.open(msaFn.c_str(), ios_base::in | ios_base::binary);
	if(!msaIn.is_open()) {
		cerr << "Unable to open MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	csfmIn.open(csfmFn.c_str(), ios_base::in | ios_base::binary);
	if(!csfmIn.is_open()) {
		cerr << "Unable to open CSFM-index '" << csfmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	hmmIn.open(hmmFn.c_str());
	if(!hmmIn.is_open()) {
		cerr << "Unable to open HMM profile '" << hmmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	ptuIn.open(ptuFn.c_str(), ios_base::in | ios_base::binary);
	if(!ptuIn.is_open()) {
		cerr << "Unable to open PTU data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* open outputs */
	if(!treeFn.empty()) {
		treeOut.open(treeFn.c_str());
		if(!treeOut.is_open()) {
			cerr << "Unable to write to tree file '" << treeFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	if(!annoFn.empty()) {
		annoOut.open(annoFn.c_str());
		if(!annoOut.is_open()) {
			cerr << "Unable to write to tree file '" << annoFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	if(!seqFn.empty()) {
		seqOut.open(seqFn.c_str());
		if(!seqOut.is_open()) {
			cerr << "Unable to write to tree file '" << seqFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		seqO.reset(&seqOut, AlphabetFactory::nuclAbc, "fasta");
	}

	/* start inspecting */
	int csLen;

	infoLog << "Inspecting MSA data ..." << endl;
	if(loadProgInfo(msaIn).bad())
		return EXIT_FAILURE;
	MSA msa;
	msa.load(msaIn);
	if(msaIn.bad()) {
		cerr << "Failed to load MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	csLen = msa.getCSLen();
	cout << "MSA loaded. Number of seq: " << msa.getNumSeq() << " CS length: " << csLen << endl;

	infoLog << "Inspecting CSFM-index ..." << endl;
	if(loadProgInfo(csfmIn).bad())
		return EXIT_FAILURE;
	CSFMIndex csfm;
	csfm.load(csfmIn);
	if(csfmIn.bad()) {
		cerr << "Failed to load CSFM-index '" << csfmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	cout << "CSFM-index loaded. Concatenated length: " << csfm.getConcatLen() << " CS length: " << csfm.getCSLen() << endl;
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
	if(loadProgInfo(ptuIn).bad())
		return EXIT_FAILURE;
	PTUnrooted ptu;
	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Unable to load Phylogenetic tree data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	const DegenAlphabet* abc = msa.getAbc();

	cout << "Phylogenetic tree loaded. Root ID: " << ptu.getRoot()->getId()
		 << " Number of leaves: " << ptu.numLeaves()
		 << " Number of nodes: " << ptu.numNodes()
		 << " Number of branches: " << ptu.numBranches()
		 << " Number of sites: " << ptu.numAlignSites() << endl;
	cout << "Overall tree log-likelihood: " << ptu.treeLoglik() << endl;

	if(showSm)
		cout << (*ptu.getModel());

	if(showDg && ptu.getDGModel() != NULL)
		cout << "Discrete Gamma Model is enabled for this tree" << endl
		     << "Number of categories used: " << ptu.getDGModel()->getK()
			 << " Shape parameter: " << ptu.getDGModel()->getShape() << endl;

	if(treeOut.is_open()) {
		infoLog << "Writing phylogenetic tree ..." << endl;
		treeOut << ptu.convertToNewickTree(nodePrefix);
	}

	if(annoOut.is_open()) {
		infoLog << "Writing tree node taxonomy annotation ..." << endl;
		for(size_t i = 0; i < ptu.numNodes(); ++i) {
			const EGriceLab::PTUnrooted::PTUNodePtr& node = ptu.getNode(i);
			annoOut << (nodePrefix + boost::lexical_cast<string>(node->getId()))
					<< "\t" << node->getTaxon() << endl;
		}
	}

	if(seqOut.is_open()) {
		infoLog << "Writing sequence alignment ..." << endl;
		for(size_t i = 0; i < ptu.numNodes(); ++i) {
			const EGriceLab::PTUnrooted::PTUNodePtr& node = ptu.getNode(i);
			if(!leafOnly || node->isLeaf())
				seqO.writeSeq(PrimarySeq(abc, nodePrefix + boost::lexical_cast<string>(node->getId()), node->getSeq().toString(), node->getTaxon()));
		}
	}

}

