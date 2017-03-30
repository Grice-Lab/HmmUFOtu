/*
 * hmmufotu-sim.cpp
 *  Generate simulated/synthesized 16S RNA reads from pre-built HmmUFOtu index
 *  Created on: Jan 10, 2017
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include <fstream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <eigen3/Eigen/Dense>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <ctime>
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_phylo.h"

using namespace std;
using namespace EGriceLab;

/**
 * default options
 */
static const string DEFAULT_FMT = "fasta";
static const string ALPHABET = "dna";

static const double DEFAULT_MEAN_LEN = 500;
static const double DEFAULT_SD_LEN = 30;
static const double DEFAULT_MIN_LEN = 0;
static const double DEFAULT_MAX_LEN = 0;

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Generate simulated multiple-aligned sequences (MSA) using a pre-built HmmUFOtu database" << endl
		 << "Usage:    " << progName << "  <HmmUFOtu-DB> <-o OUTPUT> <-N NUM-READS> [options]" << endl
		 << "Options:    -N  LONG           : number of reads to generate" << endl
	     << "            -o  FILE           : sequence OUTPUT file" << endl
		 << "            -f|--fmt  STRING   : output format [" << DEFAULT_FMT << "]" << endl
		 << "            -r|--remove-gap    : remove gaps in generated reads, seq will be unaligned" << endl
		 << "            -m|--mean-len  DBL : mean read length [" << DEFAULT_MEAN_LEN << "]" << endl
		 << "            -s|--sd-len  DBL   : standard deviation of read length [" << DEFAULT_SD_LEN << "]" << endl
		 << "            -l|--min-len  DBL  : minimum read length, 0 for no limit [" << DEFAULT_MIN_LEN << "]" << endl
		 << "            -u|--max-len  DBL  : maximum read length, 0 for no limit [" << DEFAULT_MAX_LEN << "]" << endl
		 << "            -S|--seed  INT     : random seed used for simulation, for debug purpose" << endl
		 << "            -v  FLAG           : enable verbose information" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string infn, msafn, ptufn, outfn;
	string fmt(DEFAULT_FMT);
	bool removeGap = false;
	long N = 0;
	ifstream msaIn, ptuIn;
	MSA msa;
	PTUnrooted ptu;

	double meanLen = DEFAULT_MEAN_LEN;
	double sdLen = DEFAULT_SD_LEN;
	double minLen = DEFAULT_MIN_LEN;
	double maxLen = DEFAULT_MAX_LEN;

	unsigned seed = time(NULL); // using time as default seed

	typedef boost::random::mt11213b RNG; /* random number generator type */
	typedef boost::random::uniform_int_distribution<size_t> NodeDistrib; /* node distribution in tree */
	typedef boost::random::uniform_01<> BranchDistrib; /* branching point distribution on any branch */
	typedef boost::random::uniform_smallint<> LocDistrib; /* location distribution on csLen */
	typedef boost::random::normal_distribution<> SizeDistrib; /* read size distribution */
	typedef boost::random::uniform_01<> GapDistrib; /* gap observing distribution */
	typedef boost::random::discrete_distribution<int8_t> BaseDistrib; /* base (nucleotide) distribution */
	typedef BaseDistrib::param_type BaseParam; /* base distribution parameters */

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(!(cmdOpts.numMainOpts() == 1 && cmdOpts.hasOpt("-o") && cmdOpts.hasOpt("-N"))) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	infn = cmdOpts.getMainOpt(0);
	msafn = infn + ".msa";
	ptufn = infn + ".ptu";
	msaIn.open(msafn.c_str(), ios_base::in | ios_base::binary);
	if(!msaIn.is_open()) {
		cerr << "Unable to open " << msafn << " : " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	ptuIn.open(ptufn.c_str(), ios_base::in | ios_base::binary);
	if(!ptuIn.is_open()) {
		cerr << "Unable to open " << ptufn << " : " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-N"))
		N = ::atol(cmdOpts.getOptStr("-N"));
	if(!(N > 0)) {
		cerr << "-N must be positive" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-o"))
		outfn = cmdOpts.getOpt("-o");
	else {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-r") || cmdOpts.hasOpt("--remove-gap"))
		removeGap = true;

	if(cmdOpts.hasOpt("-f"))
		fmt = StringUtils::toLower(cmdOpts.getOpt("-f"));
	if(cmdOpts.hasOpt("--fmt"))
		fmt = StringUtils::toLower(cmdOpts.getOpt("--fmt"));

	if(cmdOpts.hasOpt("-m"))
		meanLen = ::atof(cmdOpts.getOptStr("-m"));
	if(cmdOpts.hasOpt("--mean-len"))
		meanLen = ::atof(cmdOpts.getOptStr("--mean-len"));

	if(cmdOpts.hasOpt("-s"))
		sdLen = ::atof(cmdOpts.getOptStr("-s"));
	if(cmdOpts.hasOpt("--sd-len"))
		sdLen = ::atof(cmdOpts.getOptStr("--sd-len"));

	if(cmdOpts.hasOpt("-l"))
		minLen = ::atof(cmdOpts.getOptStr("-l"));
	if(cmdOpts.hasOpt("--min-len"))
		minLen = ::atof(cmdOpts.getOptStr("--min-len"));
	if(!(minLen >= 0)) {
		cerr << "-m must be non-negative" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-u"))
		maxLen = ::atof(cmdOpts.getOptStr("-u"));
	if(cmdOpts.hasOpt("--max-len"))
		maxLen = ::atof(cmdOpts.getOptStr("--max-len"));
	if(!(maxLen >= 0 && maxLen >= minLen)) {
		cerr << "-n must be non-negative and non-less than -m" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-S"))
		seed = ::atoi(cmdOpts.getOptStr("-S"));
	if(cmdOpts.hasOpt("--seed"))
		seed = ::atoi(cmdOpts.getOptStr("--seed"));
	srand(seed);

	if(cmdOpts.hasOpt("-v"))
		ENABLE_INFO();

	/* open SeqIO output */
	SeqIO seqOut(outfn, ALPHABET, fmt, SeqIO::WRITE);

	/* load input database */
	msa.load(msaIn);
	if(msaIn.bad()) {
		cerr << "Failed to load MSA data from " << msafn << endl;
		return EXIT_FAILURE;
	}
	else
		infoLog << "MSA data loaded, numSeq: " << msa.getNumSeq() << " csLen:" << msa.getCSLen() << endl;

	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Failed to load PTU data from " << ptufn << endl;
		return EXIT_FAILURE;
	}
	else
		infoLog << "Phylogenetic tree data loaded, numNode: " << ptu.numNodes() << " numSites:" << ptu.numAlignSites() << endl;

	if(msa.getCSLen() != ptu.numAlignSites()) {
		cerr << "Unmatched HmmUFOtu data files, please rebuid your database" << endl;
		return EXIT_FAILURE;
	}

	const int csLen = ptu.numAlignSites();
	const size_t numNodes = ptu.numNodes();

	/* constructor random sample generator and required distributions */
	RNG rng;
	NodeDistrib node_dist(0, numNodes - 1);
	BranchDistrib branch_dist;
	LocDistrib loc_dist(0, csLen - 1);
	SizeDistrib size_dist(meanLen, sdLen);
	GapDistrib gap_dist;
	double basePr[4] = {1, 1, 1, 1};
	Map<Vector4d> basePrMap(basePr, 4); /* use a map to access basePr indirectly */
	BaseDistrib base_dist(basePr);

	char rid[22]; // enough to hold all numbers up to 64-bits + a prefix char
	char desc[4096]; // enough to hold must descriptions
	const DegenAlphabet* abc = msa.getAbc();
	const char gapSym = abc->getGap().front();
	const PTUnrooted::ModelPtr& model = ptu.getModel();
	const Vector4d& pi = model->getPi();

	/* generating random reads */
	for(long n = 1; n <= N;) {
		/* simulate a branch */
		size_t id = node_dist(rng);
		PTUnrooted::PTUNodePtr cNode = ptu.getNode(id);
		if(cNode->isRoot()) /* no parent branch available */
			continue;
		PTUnrooted::PTUNodePtr pNode = cNode->getParent();
		double v = ptu.getBranchLength(pNode, cNode);
		double rc = branch_dist(rng);
		/* simulate a read range */
		int start = loc_dist(rng);
		double len = size_dist(rng);
		if(len < minLen)
			len = minLen;
		if(maxLen > 0 && len > maxLen)
			len = maxLen;
		int end = start + static_cast<int> (len);
		if(!(end < csLen)) /* outside consensus range */
			continue;

		/* simulate a read at [start, end] */
		sprintf(rid, "r%d", n);
		sprintf(desc, "ID=%ld->%ld;Name=\"%s\";AnnoDist=%f;Start=%d;End=%d;Len=%d;",
				cNode->getId(), pNode->getId(),
				rc < 0.5 ? cNode->getTaxa().c_str() : pNode->getTaxa().c_str(),
				rc < 0.5 ? v * rc : v * (1 - rc),
				start, end, end - start + 1);

//		PrimarySeq seq(abc, rid, "", desc);
		string seq;
		if(!removeGap)
			seq.append(start, gapSym);

		for(int j = start; j <= end; ++j) {
			bool isGap = gap_dist(rng) <= msa.gapWFrac(j);
			if(isGap && !removeGap)
				seq.push_back(gapSym);
			else {
				/* calculate the loglik of this branch point */
				Vector4d rLoglik = PTUnrooted::dot_product_scaled(model->Pr(v * rc), ptu.getBranchLoglik(cNode, pNode, j)) +
								   PTUnrooted::dot_product_scaled(model->Pr(v * (1 - rc)), ptu.getBranchLoglik(pNode, cNode, j));
				/* normalize and reset the probabilities of the base distribution */
				rLoglik.array() -= rLoglik.maxCoeff();
				basePrMap = rLoglik.array().exp();
				base_dist.param(BaseParam(basePr));
				seq.push_back(abc->decode(base_dist(rng)));
			}
		}

		if(!removeGap)
			seq.append(csLen - 1 - end, gapSym);

		seqOut.writeSeq(PrimarySeq(abc, rid, seq, desc));

		n++;
	}
}
