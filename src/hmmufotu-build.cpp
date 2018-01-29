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
 * hmmufotu-build.cpp
 * Build HmmUFOtu index files from a MSA file
 * Index files include an optional msa file, an hmm file, a csfm file and a ptu file
 *  Created on: Feb 2, 2017
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <string>
#include "HmmUFOtu.h"

#ifndef SRC_DATADIR
#define SRC_DATADIR "."
#endif

#ifndef PKG_DATADIR
#define PKG_DATADIR "."
#endif

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::HmmUFOtu;

/** default values */
static const double DEFAULT_SYMFRAC = 0.5;
static const string DEFAULT_DM_FILE = "gg_97_otus.dm";
static const string DEFAULT_SM_TYPE = "GTR";
static const string DEFAULT_SM_NAME = "gg_97_otus";
static const string ALPHABET = "dna";
static const int DEFAULT_DG_CATEGORY = 4;
static const int MIN_DG_CATEGORY = 2;
static const int MAX_DG_CATEGORY = 8;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Build an HmmUFOtu database from reference MSA and phylogenetic tree files" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <MSA-FILE> <TREE-FILE> [options]" << endl
		 << "MSA-FILE  FILE                   : multiple-sequence aligned (MSA) input" << endl
		 << "TREE-FILE  FILE                  : phylogenetic-tree file build on the MSA sequences" << endl
		 << "Options:    -n  STR              : database name (prefix), use 'MSA-FILE' by default" << endl
		 << "            --fmt  STR           : MSA format, supported format: 'fasta'" << endl
		 << "            -f|--symfrac  DOUBLE : conservation threshold for considering a site as a Match state in HMM [" << DEFAULT_SYMFRAC << "]" << endl
		 << "            -a|--anno  FILE      : use tab-delimited taxonamy annotation file for the sequences in the MSA and TREE files" << endl
		 << "            -dm  FILE            : use customized trained Dirichlet Model in FILE instead of the build-in file" << endl
		 << "            -s|--sub-model STR   : use built-in DNA Substitution Model of this type, must be one of GTR, TN93, HKY85, F81, K80 or JC69 [" << DEFAULT_SM_TYPE << "]" << endl
		 << "            -sm  FILE            : use customized trained DNA Substitution Model in FILE instead of the build-in model, will override -s if specified" << endl
		 << "            --no-hmm FLAG        : do not build the Hmm profile. Users should build the Hmm profile by 3rd party programs, i.e. HMMER3" << endl
		 << "            -V|--var FLAG        : enable among-site rate varation evaluation of the tree, using a Discrete Gamma Distribution based model" << endl
		 << "            -k INT               : number of Discrete Gamma Distribution categories to evaluate the tree, ignored if -V not set [" << DEFAULT_DG_CATEGORY << "]" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string seqFn, treeFn, dbName, annoFn;
	ifstream dmIn, smIn, seqIn, treeIn, annoIn;
	ofstream msaOut, csfmOut, hmmOut, ptuOut;
	string fmt;
	string smType = DEFAULT_SM_TYPE;
	double symfrac = DEFAULT_SYMFRAC;
	string dmFn;
	string smFn;
	bool noHmm = false;
	bool isVar = false;
	int K = DEFAULT_DG_CATEGORY;

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

	if(cmdOpts.numMainOpts() != 2) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	seqFn = cmdOpts.getMainOpt(0);
	treeFn = cmdOpts.getMainOpt(1);

	if(cmdOpts.hasOpt("-V") || cmdOpts.hasOpt("--var"))
		isVar = true;

	if(cmdOpts.hasOpt("-n"))
		dbName = cmdOpts.getOpt("-n");

	if(cmdOpts.hasOpt("--fmt"))
		fmt = cmdOpts.getOpt("--fmt");

	if(cmdOpts.hasOpt("-f"))
		symfrac = atof(cmdOpts.getOptStr("-f"));
	if(cmdOpts.hasOpt("-symfrac"))
		symfrac = atof(cmdOpts.getOptStr("-symfrac"));

	if(cmdOpts.hasOpt("-a"))
		annoFn = cmdOpts.getOpt("-a");
	if(cmdOpts.hasOpt("--anno"))
		annoFn = cmdOpts.getOpt("--anno");

	dmFn = PKG_DATADIR + string("/") + DEFAULT_DM_FILE;
	if(!ifstream(dmFn.c_str()).good())
		dmFn = SRC_DATADIR + string("/") + DEFAULT_DM_FILE;
	if(cmdOpts.hasOpt("-dm"))
		dmFn = cmdOpts.getOpt("-dm");

	if(cmdOpts.hasOpt("-s"))
		smType = cmdOpts.getOpt("-s");
	if(cmdOpts.hasOpt("--sub-model"))
		smType = cmdOpts.getOpt("--sub-model");
	smFn = PKG_DATADIR + string("/") + DEFAULT_SM_NAME + "_" + smType + SUB_MODEL_FILE_SUFFIX;
	if(!ifstream(smFn.c_str()).good())
		smFn = SRC_DATADIR + string("/") + DEFAULT_SM_NAME + "_" + smType + SUB_MODEL_FILE_SUFFIX;

	if(cmdOpts.hasOpt("-sm"))
		smFn = cmdOpts.getOpt("-sm");

	if(cmdOpts.hasOpt("--no-hmm"))
		noHmm = true;

	if(cmdOpts.hasOpt("-k"))
		K = atoi(cmdOpts.getOptStr("-k"));

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* check options */
	if(fmt.empty()) {
		if(StringUtils::endsWith(seqFn, ".fasta") || StringUtils::endsWith(seqFn, ".fas")
		|| StringUtils::endsWith(seqFn, ".fa") || StringUtils::endsWith(seqFn, ".fna"))
			fmt = "fasta";
		else {
			cerr << "Unrecognized format of MSA file '" << seqFn << "'" << endl;
			return EXIT_FAILURE;
		}
	}
	if(fmt != "fasta") {
		cerr << "Unsupported sequence format '" << fmt << "'" << endl;
		return EXIT_FAILURE;
	}

	if(!NewickTree::isNewickFileExt(treeFn)) {
		cerr << "Unrecognized TREE-FILE format, must be in Newick format" << endl;
		return EXIT_FAILURE;
	}

	if(!(MIN_DG_CATEGORY <= K && K <= MAX_DG_CATEGORY)) {
		cerr << "-k must be an integer between " << MIN_DG_CATEGORY << " and " << MAX_DG_CATEGORY << endl;
		return EXIT_FAILURE;
	}

	if(!(symfrac >= 0 && symfrac <= 1)) {
		cerr << "-f|--symfrac must between 0 and 1" << endl;
		return EXIT_FAILURE;
	}

	/* open inputs */
	seqIn.open(seqFn.c_str());
	if(!seqIn.is_open()) {
		cerr << "Unable to open '" << seqFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	dmIn.open(dmFn.c_str());
	if(!dmIn.is_open()) {
		cerr << "Unable to open '" << dmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	smIn.open(smFn.c_str());
	if(!smIn.is_open()) {
		cerr << "Unable to open '" << smFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	treeIn.open(treeFn.c_str());
	if(!treeIn.is_open()) {
		cerr << "Unable to open '" << treeFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(!annoFn.empty()) {
		annoIn.open(annoFn.c_str());
		if(!annoIn.is_open()) {
			cerr << "Unable to open '" << annoFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* set dbName */
	if(dbName.empty())
		dbName = StringUtils::basename(seqFn);
	dbName += "_" + smType + (isVar ? "_dG" : "");

	string msaFn = dbName + MSA_FILE_SUFFIX;
	string csfmFn = dbName + CSFM_FILE_SUFFIX;
	string hmmFn = dbName + HMM_FILE_SUFFIX;
	string ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;

	/* open output files */
	msaOut.open(msaFn.c_str(), ios_base::out | ios_base::binary);
	if(!msaOut.is_open()) {
		cerr << "Unable to write to '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	csfmOut.open(csfmFn.c_str(), ios_base::out | ios_base::binary);
	if(!csfmOut.is_open()) {
		cerr << "Unable to write to '" << csfmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	hmmOut.open(hmmFn.c_str());
	if(!hmmOut.is_open()) {
		cerr << "Unable to write to '" << hmmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	ptuOut.open(ptuFn.c_str(), ios_base::out | ios_base::binary);
	if(!ptuOut.is_open()) {
		cerr << "Unable to write to '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* build msa */
	MSA msa;
	if(msa.loadMSA(ALPHABET, seqIn, fmt) >= 0)
		infoLog << "MSA loaded" << endl;
	else {
		cerr << "Unable to load MSA from '" << seqFn << "'" << endl;
		return EXIT_FAILURE;
	}
	msa.setName(dbName);
	msa.prune();
	infoLog << "MSA pruned" << endl;
	infoLog << "MSA database created for " << msa.getNumSeq() << " X " << msa.getCSLen() << " aligned sequences" << endl;

	/* build csfm */
	CSFMIndex csfm;
	csfm.build(msa);
	if(csfm.isInitiated())
		infoLog << "CSFM index built" << endl;
	else {
		cerr << "Unable to build CSFM index" << endl;
		return EXIT_FAILURE;
	}

	/* build hmm, if requested */
	/* Load in BandedHmmPrior for the HMM training */
	BandedHMMP7Prior hmmPrior;
	dmIn >> hmmPrior;
	if(dmIn.bad()) {
		cerr << "Failed to read in the HMM Prior file '" << dmFn << "'" << endl;
		return EXIT_FAILURE;
	}
	BandedHMMP7 hmm; /* construct an empty profile */
	hmm.setName(dbName);
	hmm.setHmmVersion(getProgFullName(progName, progVer));
	if(!noHmm) {
		hmm.build(msa, symfrac, hmmPrior);
		infoLog << "Banded HMM profile trained" << endl;
	}

	/* build ptu */
	NewickTree NTree;
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
		cerr << "Unmatched MSA and Tree. Found " << nRead << " leaf sequences from MSA but expecting " << nLeaves << " leaves in the Phylogenetic Tree " << endl;
		return EXIT_FAILURE;
	}
	else
		infoLog << "MSA loaded into Phylogenetic Tree" << endl;
	if(annoIn.is_open()) {
		tree.loadAnnotation(annoIn);
		if(annoIn.bad()) {
			cerr << "Failed to load taxonomy annotation from '" << annoFn << "'" << endl;
			return EXIT_FAILURE;
		}
		else
			infoLog << "Taxonomy annotation loaded" << endl;
	}
	tree.formatName();
	infoLog << "Taxon names formatted" << endl;

	tree.annotate();
	infoLog << "Unnamed tree nodes annotated" << endl;
//	for(int i = 0; i < tree.numNodes(); ++i)
//		cerr << "ID: " << tree.getNode(i)->getId() << " annotation: " << tree.getNode(i)->getTaxon() << endl;

	tree.calcNodeHeight();
	infoLog << "Node height calculated" << endl;

	/* read DNA sub model from model file */
	DNASubModel* model = NULL;
	string line, tag, type;
	while(smIn >> tag) {
		if(tag[0] == '#') { /* comment or header */
			std::getline(smIn, line); /* ignore the entire line */
			continue;
		}
		if(tag == "Type:") {
			smIn >> type; // read in model type
			model = DNASubModelFactory::createModel(type); /* dynamic construction */
			smIn >> *model; /* read the remaining file */
		}
	}

	if(model != NULL)
		infoLog << "DNA Substitution Model loaded" << endl;
	else {
		cerr << "Unable to load DNA Substitution Model from: '" << smFn << "'" << endl;
		return EXIT_FAILURE;
	}
	tree.setModel(model);

	/* initiation the tree costs */
	tree.initRootLoglik();
	tree.initBranchLoglik();
	tree.initLeafMat();

	/* make initial evaluation at the original root */
	const PTUnrooted::PTUNodePtr root = tree.getRoot();
	if(!isVar)
		infoLog << "Evaluating Phylogenetic Tree at root id: " << root->getId() << endl;
	else
		infoLog << "Evaluating Phylogenetic Tree at root id: " << root->getId() << " with fixed rate model first" << endl;
	tree.evaluate(); /* only evaluate, do not cache the root loglik */
//	infoLog << "tree log-liklihood: " << tree.treeLoglik() << endl;

	/* construct DG model, if isVar is set */
	if(isVar) {
		infoLog << "Estimating the shape parameter of the Discrete Gamma Distributin based among-site variation ..." << endl;
		VectorXd numMut(tree.numAlignSites());
		for(int j = 0; j < tree.numAlignSites(); ++j)
			numMut(j) = tree.estimateNumMutations(j);
		double alpha = DiscreteGammaModel::estimateShape(numMut);
		if(alpha == inf)
			cerr << "Unable to estimate the shape parameter with less than 2 alignment sites" << endl;
		else if(alpha <= 0)
			cerr << "Unable to estimate the shape parameter with near invariant rates, reducing to fixed rate model" << endl;
		else {
			infoLog << "Estimated alpha = " << alpha << endl;
			tree.setDGModel(DiscreteGammaModel(K, alpha));
			tree.resetBranchLoglik(); /* reset all cached values */
			tree.resetRootLoglik();   /* reset cached root value */
		}
	}

	if(!isVar)
		infoLog << "Evaluating Phylogenetic Tree at all other " << (tree.numNodes() - 1) << " nodes" << endl;
	else
		infoLog << "Re-evaluating Phylogenetic Tree at all " << tree.numNodes() << " nodes" << endl;

	const size_t numNodes = tree.numNodes();
	for(size_t i = 0; i < numNodes; ++i) {
//		cerr << "Setting root at " << i << endl;
		tree.setRoot(i);
		tree.evaluate();
//		cerr << "loglik: " << tree.loglik() << endl;
	}
	/* reset to original root, and evaluate its root Loglik */
	tree.setRoot(root);
	tree.loglik(); /* evaluate the root and store value */
	infoLog << "Final Tree log-liklihood: " << tree.treeLoglik() << endl;

	/* infer the ancestor seq of all intermediate nodes */
	tree.inferSeq();
	infoLog << "Ancestor sequence of all intermediate nodes inferred" << endl;

	infoLog << "Saving database files ..." << endl;
	/* write database files, all with prepend program info */
	saveProgInfo(msaOut);
	msa.save(msaOut);
	if(msaOut.bad()) {
		cerr << "Unable to save MSA: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "MSA saved" << endl;

	saveProgInfo(csfmOut);
	csfm.save(csfmOut);
	if(csfmOut.bad()) {
		cerr << "Unable to save CSFM index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "CSFM saved" << endl;

	if(!noHmm) {
		hmmOut << hmm;
		if(hmmOut.bad()) {
			cerr << "Unable to save HMM profile: " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		infoLog << "Banded HMM profile saved" << endl;
	}

	saveProgInfo(ptuOut);
	tree.save(ptuOut);
	if(ptuOut.bad()) {
		cerr << "Unable to save Phylogenetic Tree index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Phylogenetic Tree index saved" << endl;
}
