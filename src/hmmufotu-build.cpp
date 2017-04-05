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

#ifndef DM_DATADIR
#define DM_DATADIR "."
#endif

#ifndef SM_DATADIR
#define SM_DATADIR "."
#endif

using namespace std;
using namespace EGriceLab;

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
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Build a HmmUFOtu database" << endl
		 << "Usage:    " << progName << "  <MSA-FILE> <TREE-FILE> [options]" << endl
		 << "MSA-FILE  FILE                   : multiple-sequence aligned (MSA) input" << endl
		 << "TREE-FILE  FILE                  : phylogenetic-tree file build on the MSA sequences" << endl
		 << "Options:    -n  STR              : database name (prefix), use 'MSA-FILE' by default" << endl
		 << "            -f|--symfrac  DOUBLE : conservation threshold for considering a site as a Match state in HMM [" << DEFAULT_SYMFRAC << "]" << endl
		 << "            -a|--anno  FILE      : use tab-delimited taxonamy annotation file for the sequences in the MSA and TREE files" << endl
		 << "            -dm  FILE            : use customized trained Dirichlet Model in FILE instead of the build-in file" << endl
		 << "            -s|--sub-model STR   : use built-in DNA Substitution Model of this type, must be one of GTR, TN93, HKY85, F81, K80 or JC69 [" << DEFAULT_SM_TYPE << "]" << endl
		 << "            -sm  FILE            : use customized trained DNA Substitution Model in FILE instead of the build-in model, will override -s if specified" << endl
		 << "            --no-hmm FLAG        : do not build the Hmm profile. Users should build the Hmm profile by 3rd party programs, i.e. HMMER3" << endl
		 << "            -V|--var FLAG        : enable among-site rate varation evaluation of the tree, using a Discrete Gamma Distribution based model" << endl
		 << "            -k INT               : number of Discrete Gamma Distribution categories to evaluate the tree, ignored if -V not set [" << DEFAULT_DG_CATEGORY << "]" << endl
		 << "            -v  FLAG             : enable verbose information" << endl
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
	string dmFn = DM_DATADIR + string("/") + DEFAULT_DM_FILE;
	string smFn;
	bool noHmm = false;
	bool isVar = false;
	int K = DEFAULT_DG_CATEGORY;

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

	seqFn = cmdOpts.getMainOpt(0);
	treeFn = cmdOpts.getMainOpt(1);

	if(cmdOpts.hasOpt("-V") || cmdOpts.hasOpt("--var"))
		isVar = true;

	if(cmdOpts.hasOpt("-n"))
		dbName = cmdOpts.getOpt("-n");

	if(cmdOpts.hasOpt("-f"))
		symfrac = atof(cmdOpts.getOptStr("-f"));
	if(cmdOpts.hasOpt("-symfrac"))
		symfrac = atof(cmdOpts.getOptStr("-symfrac"));

	if(cmdOpts.hasOpt("-a"))
		annoFn = cmdOpts.getOpt("-a");
	if(cmdOpts.hasOpt("--anno"))
		annoFn = cmdOpts.getOpt("--anno");

	if(cmdOpts.hasOpt("-dm"))
		dmFn = cmdOpts.getOpt("-dm");

	if(cmdOpts.hasOpt("-s"))
		smType = cmdOpts.getOpt("-s");
	if(cmdOpts.hasOpt("--sub-model"))
		smType = cmdOpts.getOpt("--sub-model");
	smFn = SM_DATADIR + string("/") + DEFAULT_SM_NAME + "_" + smType + SUB_MODEL_FILE_SUFFIX;

	if(cmdOpts.hasOpt("-sm"))
		smFn = cmdOpts.getOpt("-sm");

	if(cmdOpts.hasOpt("--no-hmm"))
		noHmm = true;

	if(cmdOpts.hasOpt("-k"))
		K = atoi(cmdOpts.getOptStr("-k"));

	if(cmdOpts.hasOpt("-v"))
		ENABLE_INFO();

	/* check options */
	if(StringUtils::endsWith(seqFn, ".fasta") || StringUtils::endsWith(seqFn, ".fas")
		|| StringUtils::endsWith(seqFn, ".fa") || StringUtils::endsWith(seqFn, ".fna"))
		fmt = "fasta";
	else {
		cerr << "Unrecognized format of MSA file '" << seqFn << "'" << endl;
		return EXIT_FAILURE;
	}

	if(!(StringUtils::endsWith(treeFn, ".tree"))) {
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
	if(msa.loadMSAFile(ALPHABET, seqFn, fmt) >= 0)
		infoLog << "MSA loaded" << endl;
	else {
		cerr << "Unable to load MSA from '" << seqFn << "'" << endl;
		return EXIT_FAILURE;
	}
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
	if(!noHmm) {
		hmm.build(msa, symfrac, hmmPrior);
		infoLog << "Banded HMM profile trained" << endl;
	}

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
	infoLog << "Taxa names formatted" << endl;

	tree.annotate();
	infoLog << "Unnamed tree nodes annotated" << endl;
//	for(int i = 0; i < tree.numNodes(); ++i)
//		cerr << "ID: " << tree.getNode(i)->getId() << " annotation: " << tree.getNode(i)->getTaxa() << endl;

	/* read DNA sub model from model file */
	DNASubModel* model = NULL;
	string line, tag, type;
	while(smIn >> tag) {
		if(tag.front() == '#') { /* comment or header */
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
	const PTUnrooted::PTUNodePtr& root = tree.getRoot();
	infoLog << "Evaluating Phylogenetic Tree at root id: " << root->getId() << endl;
	tree.loglik(); /* use loglik() instead of evaluate() to cache root loglik */
//	infoLog << "tree log-liklihood: " << tree.treeLoglik() << endl;

	/* estimate the shape parameter, if using DG model */
	if(isVar) {
		infoLog << "Estimating the shape parameter of the Discrete Gamma Distributin based among-site variation ..." << endl;
		VectorXi numMut(tree.numAlignSites());
		for(int j = 0; j < tree.numAlignSites(); ++j)
			numMut(j) = tree.estimateNumMutations(j);
		double alpha = DiscreteGammaModel::estimateShape(numMut);
		if(alpha == inf)
			cerr << "Unable to estimate the shape parameter with less than 2 alignment sites" << endl;
		else if(alpha < 0)
			cerr << "Unable to estimate the shape parameter with near invariant rates, reducing to fixed rate model" << endl;
		else {
			infoLog << "Estimated alpha = " << alpha << endl;
			tree.setDGModel(DiscreteGammaModel(K, alpha));
		}
	}

	infoLog << (!isVar ? "Evaluating Phylogenetic Tree on all other nodes ..."
			: "Re-evaluating Phylogenetic Tree rooting at all possible nodes ...") << endl;
	EGriceLab::PTUnrooted::PTUNodePtr oldRoot = tree.getRoot();
	size_t numNodes = tree.numNodes();
	for(size_t i = 0; i < numNodes; ++i) {
		tree.setRoot(i);
		tree.evaluate();
//		infoLog << (i + 1) << "/" << numNodes << " nodes evaluated\r";
	}
	/* reset root to original */
	tree.setRoot(oldRoot);
//	infoLog << endl << "Tree log-liklihood: " << tree.treeLoglik() << endl;
	infoLog << "Saving database files ..." << endl;

	/* write database files */
	if(!msa.save(msaOut)) {
		cerr << "Unable to save MSA: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "MSA saved" << endl;

	if(!csfm.save(csfmOut)) {
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

	if(!tree.save(ptuOut)) {
		cerr << "Unable to save Phylogenetic Tree index: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Phylogenetic Tree index saved" << endl;
}
