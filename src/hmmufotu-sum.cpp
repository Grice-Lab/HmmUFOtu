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
static const string ALIGN_FORMAT = "fasta";
static const double DEFAULT_EFFN = 10;
static const Eigen::IOFormat otuFmt(StreamPrecision, 0, "\t");

/** struct to store observed count information for a node */
struct NodeObsData {
	NodeObsData(int L, int S) :
		freq(4, L), count(S)
	{
		freq.setZero();
		count.setZero();
	}

	Matrix4Xd freq; /* observed aggregate base frequency over all samples */
	VectorXi count; /* observed sequence count for each sample separately */
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
	cerr << "Usage:    " << progName << "  <HmmUFOtu-DB> <PLACE-FILE> [PLACE-FILE2 ...]> <-o OTU-OUT> <-c CS-OUT> [options]" << endl
		 << "PLACE-FILE      FILE           : placement file(s) from hmmufotu" << endl
		 << "Options:    -o  FILE           : OTU summary output" << endl
		 << "            -c  FILE           : OTU Consensus Sequence (CS) output" << endl
		 << "            -e|--effN  double  : effective number of sequences underlying the pre-built Phylogenetic Tree used for Bayesian inference of the OTU CS with Dirichelet Density models [" << DEFAULT_EFFN << "]" << endl
		 << "            -v  FLAG           : enable verbose information" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, msaFn, ptuFn;
	vector<string> inFiles;
	vector<string> sampleNames;
	string otuFn, csFn;
	ifstream msaIn, ptuIn;
	ofstream otuOut;
	SeqIO csOut;

	double effN = DEFAULT_EFFN;

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
	else {
		cerr << "-c must be specified" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-e"))
		effN = ::atof(cmdOpts.getOptStr("-e"));
	if(cmdOpts.hasOpt("--effN"))
		effN = ::atof(cmdOpts.getOptStr("--effN"));

	if(cmdOpts.hasOpt("-v"))
		ENABLE_INFO();

	/* validate options */
	if(!(effN > 0)) {
		cerr << "-e|--effN must be positive" << endl;
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

	csOut.open(csFn, "dna", "fasta", SeqIO::WRITE);
	if(!csOut.is_open()) {
		cerr << "Unable to write to '" << csFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	/* loading database files */
	MSA msa;
	msa.load(msaIn);
	if(msaIn.bad()) {
		cerr << "Failed to load MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	int csLen = msa.getCSLen();
	infoLog << "MSA loaded." << endl;

	PTUnrooted ptu;
	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Unable to load Phylogenetic tree data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Phylogenetic tree loaded." << endl;

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
			iss.ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t');
			iss >> taxa_id;
			iss.ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t').ignore(maxIgnore, '\t');
			iss >> aln;
			if(taxa_id < 0)
				continue; /* invalid placement */
			const PTUnrooted::PTUNodePtr& node = ptu.getNode(taxa_id);
			if(otuData.find(node) == otuData.end()) /* not initiated */
				otuData.insert(std::make_pair(node, NodeObsData(L, S)));
			NodeObsData& data = otuData.find(node)->second;
			data.count(s)++;
			for(int j = 0; j < L; ++j) {
				int8_t b = abc->encode(aln[j]);
				if(b >= 0)
					data.freq(b, s)++;
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
		NodeObsData& data = otuData.find(node)->second;
		otuOut << node->getId() << "\t";
		for(int s = 0; s < S; ++s)
			otuOut << data.count(s) << "\t";
		otuOut << node->getTaxa() << endl;
		DigitalSeq seq = ptu.inferPostCS(node, data.freq, effN);
		string desc = "DBName="
				+ dbName + ";Taxonomy=\"" + node->getTaxa()
				+ "\";AnnoDist=" + boost::lexical_cast<string>(node->getAnnoDist())
				+ ";ReadCount=" + boost::lexical_cast<string>(data.count.sum())
				+ ";SampleHits=" + boost::lexical_cast<string>((data.count.array() > 0).count());
		csOut.writeSeq(PrimarySeq(seq.getAbc(), seq.getName(), seq.toString(), desc));
	}
}
