/*
 * hmmufotu-sum.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: Qi Zheng
 *      Version: v1.1
 *
 */

#include <iostream>
#include <fstream>
#include <strstream>
#include <cctype>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <limits>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp> /* for boost string split */
#include <boost/lexical_cast.hpp>
#include "HmmUFOtu.h"

using namespace std;
using namespace EGriceLab;
using namespace Eigen;

/* default values */
static const size_t maxIgnore = numeric_limits<streamsize>::max();
static const string ALPHABET = "dna";
static const string ALIGN_FORMAT = "fasta";
static const double DEFAULT_EFFN = 2;
static const int DEFAULT_MIN_NREAD = 1;
static const int DEFAULT_MIN_NSAMPLE = 1;
//static const Eigen::IOFormat otuFmt(StreamPrecision, 0, "\t");

/** struct to store observed count information for a node */
struct NodeObsData {
	NodeObsData(int L, int S) :
		freq(4, L), gap(L), count(S)
	{
		freq.setZero();
		gap.setZero();
		count.setZero();
	}

	Matrix4Xd freq;  /* observed aggregate base frequency over all samples */
	RowVectorXd gap; /* observed aggregate gap over all samples */
	VectorXi count;  /* observed sequence count for each sample separately */
};

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Get OTU summary table and taxonomy informatin based on hmmufotu placement and alignment files" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <HmmUFOtu-DB> <PLACE-FILE> [PLACE-FILE2 ...]> <-o OTU-OUT> [options]" << endl
		 << "PLACE-FILE      FILE           : placement file(s) from hmmufotu" << endl
		 << "Options:    -o  FILE           : OTU summary output" << endl
		 << "            -c  FILE           : OTU Consensus Sequence (CS) alignment output" << endl
		 << "            -t  FILE           : OTU tree output" << endl
		 << "            -e|--effN  double  : effective number of sequences (pseudo-count) underlying the pre-built Phylogenetic Tree used for Bayesian inferencing OTUs' CS with Dirichelet Density models [" << DEFAULT_EFFN << "]" << endl
		 << "            -n  INT      : minimum number of observed reads required to define an OTU across all samples [" << DEFAULT_MIN_NREAD << "]" << endl
		 << "            -s  INT   : minimum number of observed samples required to define an OTU" << DEFAULT_MIN_NSAMPLE << "]" << endl
		 << "            -v  FLAG           : enable verbose information" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, msaFn, ptuFn;
	vector<string> inFiles;
	vector<string> sampleNames;
	string otuFn, csFn, treeFn;
	ifstream msaIn, ptuIn;
	ofstream otuOut, treeOut;
	SeqIO csOut;

	double effN = DEFAULT_EFFN;
	int minRead = DEFAULT_MIN_NREAD;
	int minSample = DEFAULT_MIN_NSAMPLE;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(!(cmdOpts.numMainOpts() > 1)) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	dbName = cmdOpts.getMainOpt(0);
	for(int i = 1; i < cmdOpts.numMainOpts(); ++i) {
		inFiles.push_back(cmdOpts.getMainOpt(i));
		sampleNames.push_back(StringUtils::basename(inFiles.back()));
	}

	if(cmdOpts.hasOpt("-o"))
		otuFn = cmdOpts.getOpt("-o");
	else {
		cerr << "-o must be specified" << endl;
		return EXIT_FAILURE;
	}
	if(cmdOpts.hasOpt("-c"))
		csFn = cmdOpts.getOpt("-c");
	if(cmdOpts.hasOpt("-t"))
		treeFn = cmdOpts.getOpt("-t");

	if(cmdOpts.hasOpt("-e"))
		effN = ::atof(cmdOpts.getOptStr("-e"));
	if(cmdOpts.hasOpt("--effN"))
		effN = ::atof(cmdOpts.getOptStr("--effN"));

	if(cmdOpts.hasOpt("-n"))
		minRead = ::atoi(cmdOpts.getOptStr("-n"));
	if(cmdOpts.hasOpt("-s"))
		minSample = ::atoi(cmdOpts.getOptStr("-s"));

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* validate options */
	if(!(effN > 0)) {
		cerr << "-e|--effN must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(!(minRead > 0)) {
		cerr << "-n must be positive integer" << endl;
		return EXIT_FAILURE;
	}
	if(!(minSample > 0)) {
		cerr << "-s must be positive integer" << endl;
		return EXIT_FAILURE;
	}

	/* set filenames */
	msaFn = dbName + MSA_FILE_SUFFIX;
	ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;

	/* open inputs */
	msaIn.open(msaFn.c_str(), ios_base::in | ios_base::binary);
	if(!msaIn) {
		cerr << "Unable to open MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	ptuIn.open(ptuFn.c_str(), ios_base::in | ios_base::binary);
	if(!ptuIn) {
		cerr << "Unable to open PTU data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* open outputs */
	otuOut.open(otuFn.c_str());
	if(!otuOut.is_open()) {
		cerr << "Unable to write to '" << otuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(!csFn.empty()) {
		csOut.open(csFn, ALPHABET, ALIGN_FORMAT, SeqIO::WRITE);
		if(!csOut.is_open()) {
			cerr << "Unable to write to '" << csFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	if(!treeFn.empty()) {
		treeOut.open(treeFn.c_str());
		if(!treeOut.is_open()) {
			cerr << "Unable to write to '" << treeFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* loading database files */
	MSA msa;
	msa.load(msaIn);
	if(msaIn.bad()) {
		cerr << "Failed to load MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	int csLen = msa.getCSLen();
	infoLog << "MSA loaded" << endl;

	PTUnrooted ptu;
	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Unable to load Phylogenetic tree data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Phylogenetic tree loaded" << endl;
	ptu.setRoot(0);

	const DegenAlphabet* abc = msa.getAbc();
	const int S = inFiles.size();
	const int L = ptu.numAlignSites();
	const size_t N = ptu.numNodes();

	infoLog << "Getting sample names from placement file names" << endl;

	boost::unordered_map<PTUnrooted::PTUNodePtr, NodeObsData> otuData;

	for(int s = 0; s < S; ++s) {
		infoLog << "Processing " << inFiles[s] << " ..." << endl;
		ifstream in(inFiles[s].c_str());
		if(!in.is_open()) {
			cerr << "Unable to open '" << inFiles[s] << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		string line;
		long taxa_id;
		string aln;
		while((std::getline(in, line))) {
			if(StringUtils::startsWith(line, "id"))
				continue;
			istringstream iss(line);
			iss.ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t');
			iss >> aln;
			iss.ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t');
			iss >> taxa_id;

			if(taxa_id < 0)
				continue; /* invalid placement */
			const PTUnrooted::PTUNodePtr& node = ptu.getNode(taxa_id);
			if(otuData.find(node) == otuData.end()) /* not initiated */
				otuData.insert(std::make_pair(node, NodeObsData(L, S)));
			NodeObsData& data = otuData.find(node)->second;
			data.count(s)++;
			for(int j = 0; j < L; ++j) {
				int8_t b = abc->encode(::toupper(aln[j]));
				if(b >= 0)
					data.freq(b, j)++;
				else
					data.gap(j)++;
			}
		}
	}

	/* output OTU table and alignment */
	infoLog << "Generating output files" << endl;
	otuOut << "OTU_id\t" << boost::join(sampleNames, "\t") << "\tTaxonomy" << endl;
	for(int i = 0; i < N; ++i) {
		const PTUnrooted::PTUNodePtr& node = ptu.getNode(i);
		if(otuData.find(node) == otuData.end())
			continue;
		otuData.begin();
		NodeObsData& data = otuData.find(node)->second;
		int nRead = data.count.sum();
		int nSample = (data.count.array() > 0).count();
		if(!(nRead >= minRead && nSample >= minSample)) /* filter OTUs */
			continue;

		otuOut << node->getId() << "\t";
		for(int s = 0; s < S; ++s)
			otuOut << data.count(s) << "\t";
		otuOut << node->getTaxa() << endl;

		if(csOut.is_open()) {
			DigitalSeq seq = ptu.inferPostCS(node, data.freq, data.gap, effN);
			string desc = "DBName="
					+ dbName + ";Taxonomy=\"" + node->getTaxa()
					+ "\";AnnoDist=" + boost::lexical_cast<string>(node->getAnnoDist())
					+ ";ReadCount=" + boost::lexical_cast<string>(nRead)
					+ ";SampleHits=" + boost::lexical_cast<string>(nSample);
			csOut.writeSeq(PrimarySeq(seq.getAbc(), seq.getName(), seq.toString(), desc));
		}
	}

	/* write the tree */
	if(treeOut.is_open())
		ptu.exportTree(treeOut, otuData, "newick");
}
