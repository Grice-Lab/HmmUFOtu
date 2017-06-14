/*
 * hmmufotu-sim.cpp
 *  Generate simulated/synthesized 16S RNA reads from pre-built HmmUFOtu index
 *  Created on: Jan 10, 2017
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/algorithm/string.hpp> /* for boost string split */
#include <Eigen/Dense>
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

static const double DEFAULT_MAX_DIST = EGriceLab::inf;
static const double DEFAULT_MEAN_SIZE = 500;
static const double DEFAULT_SD_SIZE = 30;
static const double DEFAULT_MIN_SIZE = 0;
static const double DEFAULT_MAX_SIZE = 0;
static const int DEFAULT_READ_LEN = -1;
static const size_t MAX_IGNORE = numeric_limits<streamsize>::max();
static const char GAP_SYM = '-';
static const char PAD_SYM = '.';

/* a simple POD for holding a CS locus */
struct CSLoc {
	/* constructors */
	CSLoc() : csStart(), csEnd() { }
	CSLoc(int csStart, int csEnd) : csStart(csStart), csEnd(csEnd) { }

	int csStart; /* 1-based CS start */
	int csEnd;   /* 1-based CS end */
};

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Generate simulated multiple-aligned sequences (MSA) using a pre-built HmmUFOtu database" << endl
		 << "Usage:    " << progName << "  <HmmUFOtu-DB> <SEQ-OUT> [MATE-OUT] <-N NUM-READS> [options]" << endl
		 << "Options:    SEQ-OUT  FILE       : OUTPUT file" << endl
		 << "            MATE-OUT  FILE      : optional OUTPUT file for paired-end mode, this will suppress -k|--keep-gap option" << endl
		 << "            -N  LONG            : number of reads/pairs to generate" << endl
		 << "            -f|--fmt  STRING    : output format [" << DEFAULT_FMT << "]" << endl
		 << "            -k|--keep-gap FLAG  : keep simulated gaps in generated reads, so final seq will be aligned" << endl
		 << "            -d|--max-dist       : maximum height allowed for simulated reads (as shorted phylogenetic distance to any leaf) [" << DEFAULT_MAX_DIST << "]" << endl
		 << "            -m|--mean-size  DBL : mean 16S amplicon size [" << DEFAULT_MEAN_SIZE << "]" << endl
		 << "            -s|--sd-size  DBL   : standard deviation of 16S amplicon size [" << DEFAULT_SD_SIZE << "]" << endl
		 << "            -l|--min-size  DBL  : minimum 16S amplicon size, 0 for no limit [" << DEFAULT_MIN_SIZE << "]" << endl
		 << "            -u|--max-size  DBL  : maximum 16S amplicon size, 0 for no limit [" << DEFAULT_MAX_SIZE << "]" << endl
		 << "            -r|--read-len  INT  : read length for generating single/paired-end reads, set to -1 to use the actual amplicon size [" << DEFAULT_READ_LEN << "]" << endl
		 << "            -R|--region  STRING : BED file for restricted consensus region where simulated reads should be drawn; setting this will ignore -m,-s,-l,-u togather" << endl
		 << "            -S|--seed  INT      : random seed used for simulation, for debug purpose" << endl
		 << "            -v  FLAG            : enable verbose information" << endl
		 << "            -h|--help           : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string inFn, msaFn, ptuFn, outFn, mateFn, regionFn;
	bool keepGap = false;
	long N = 0;
	ifstream msaIn, ptuIn, regionIn;
	MSA msa;
	PTUnrooted ptu;

	double maxDist = DEFAULT_MAX_DIST;
	double meanSize = DEFAULT_MEAN_SIZE;
	double sdSize = DEFAULT_SD_SIZE;
	double minSize = DEFAULT_MIN_SIZE;
	double maxSize = DEFAULT_MAX_SIZE;
	int readLen = DEFAULT_READ_LEN;
	vector<CSLoc> myLoci;

	unsigned seed = time(NULL); // using time as default seed

	typedef boost::random::mt11213b RNG; /* random number generator type */
	typedef boost::random::discrete_distribution<size_t> NodeDistrib; /* node distribution in tree */
	typedef boost::random::uniform_01<> BranchDistrib; /* branching point distribution on any branch */
	typedef boost::random::uniform_smallint<> LocDistrib; /* location distribution either on csLen, or myLoci, if provided */
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

	if(!((cmdOpts.numMainOpts() == 2 || cmdOpts.numMainOpts() == 3) && cmdOpts.hasOpt("-N"))) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	inFn = cmdOpts.getMainOpt(0);

	outFn = cmdOpts.getMainOpt(1);
	if(cmdOpts.numMainOpts() == 3)
		mateFn = cmdOpts.getMainOpt(2);

	if(cmdOpts.hasOpt("-N"))
		N = ::atol(cmdOpts.getOptStr("-N"));

	if(cmdOpts.hasOpt("-k") || cmdOpts.hasOpt("--keep-gap"))
		keepGap = true;

	if(cmdOpts.hasOpt("-d"))
		maxDist = ::atof(cmdOpts.getOptStr("-d"));
	if(cmdOpts.hasOpt("--max-dist"))
		maxDist = ::atof(cmdOpts.getOptStr("--max-dist"));

	if(cmdOpts.hasOpt("-m"))
		meanSize = ::atof(cmdOpts.getOptStr("-m"));
	if(cmdOpts.hasOpt("--mean-size"))
		meanSize = ::atof(cmdOpts.getOptStr("--mean-size"));

	if(cmdOpts.hasOpt("-s"))
		sdSize = ::atof(cmdOpts.getOptStr("-s"));
	if(cmdOpts.hasOpt("--sd-len"))
		sdSize = ::atof(cmdOpts.getOptStr("--sd-size"));

	if(cmdOpts.hasOpt("-l"))
		minSize = ::atof(cmdOpts.getOptStr("-l"));
	if(cmdOpts.hasOpt("--min-size"))
		minSize = ::atof(cmdOpts.getOptStr("--min-size"));

	if(cmdOpts.hasOpt("-u"))
		maxSize = ::atof(cmdOpts.getOptStr("-u"));
	if(cmdOpts.hasOpt("--max-size"))
		maxSize = ::atof(cmdOpts.getOptStr("--max-size"));

	if(cmdOpts.hasOpt("-r"))
		readLen = ::atoi(cmdOpts.getOptStr("-r"));
	if(cmdOpts.hasOpt("--read-len"))
		readLen = ::atoi(cmdOpts.getOptStr("--read-len"));

	if(cmdOpts.hasOpt("-R"))
		regionFn = cmdOpts.getOpt("-R");
	if(cmdOpts.hasOpt("--region"))
		regionFn = cmdOpts.getOpt("--region");

	if(cmdOpts.hasOpt("-S"))
		seed = ::atoi(cmdOpts.getOptStr("-S"));
	if(cmdOpts.hasOpt("--seed"))
		seed = ::atoi(cmdOpts.getOptStr("--seed"));

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* validate options */
	if(!(N > 0)) {
		cerr << "-N must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(!(meanSize > 0)) {
		cerr << "-m|--min-size must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(!(sdSize > 0)) {
		cerr << "-s|--sd-size must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(!(minSize >= 0)) {
		cerr << "-l|--min-size must be non-negative" << endl;
		return EXIT_FAILURE;
	}
	if(!(maxSize >= 0 && maxSize >= minSize)) {
		cerr << "-u|--max-size must be non-negative and non-less than -l|--min-size" << endl;
		return EXIT_FAILURE;
	}

	/* open inputs */
	msaFn = inFn + ".msa";
	ptuFn = inFn + ".ptu";
	msaIn.open(msaFn.c_str(), ios_base::in | ios_base::binary);
	if(!msaIn.is_open()) {
		cerr << "Unable to open " << msaFn << " : " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	ptuIn.open(ptuFn.c_str(), ios_base::in | ios_base::binary);
	if(!ptuIn.is_open()) {
		cerr << "Unable to open " << ptuFn << " : " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(!regionFn.empty()) {
		regionIn.open(regionFn.c_str(), ios_base::in);
		if(!regionIn.is_open()) {
			cerr << "Unable to open " << regionFn << " : " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* open outputs */
	SeqIO seqOut(outFn, AlphabetFactory::nuclAbc, DEFAULT_FMT, SeqIO::WRITE, -1);
	if(!seqOut.is_open()) {
		cerr << "Unable to write seq to '" << outFn << "' : " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	SeqIO mateOut;
	if(!mateFn.empty()) {
		keepGap = false; /* suppress -k if paired end */
		mateOut.open(mateFn, AlphabetFactory::nuclAbc, DEFAULT_FMT, SeqIO::WRITE, -1);
		if(!mateOut.is_open()) {
			cerr << "Unable to write mate to '" << mateFn << "' : " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* load input database */
	msa.load(msaIn);
	if(msaIn.bad()) {
		cerr << "Failed to load MSA data from " << msaFn << endl;
		return EXIT_FAILURE;
	}
	else
		infoLog << "MSA data loaded, numSeq: " << msa.getNumSeq() << " csLen:" << msa.getCSLen() << endl;

	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Failed to load PTU data from " << ptuFn << endl;
		return EXIT_FAILURE;
	}
	else
		infoLog << "Phylogenetic tree data loaded, numNode: " << ptu.numNodes() << " numSites:" << ptu.numAlignSites() << endl;

	if(msa.getCSLen() != ptu.numAlignSites()) {
		cerr << "Unmatched HmmUFOtu data files, please rebuild your database" << endl;
		return EXIT_FAILURE;
	}

	const int csLen = ptu.numAlignSites();
	const size_t numNodes = ptu.numNodes();

	/* read restricted regions, if provided */
	if(regionIn.is_open()) {
		int start, end;
		while(regionIn) {
			regionIn.ignore(MAX_IGNORE, '\t') >> start;
			regionIn.ignore(MAX_IGNORE, '\t') >> end;
			regionIn.ignore(MAX_IGNORE, '\n');
			if(!(0 <= start && start < end && end <= csLen)) {
				warningLog << "provided region (" << start << "," << end << "] is not in the consensus sequence range, ignored" << endl;
				continue;
			}
			myLoci.push_back(CSLoc(start + 1, end));
			if(regionIn.peek() == EOF)
				break;
		}
		infoLog << "Read in " << myLoci.size() << " restricted regions" << endl;
	}

	/* constructor random sample generator and required distributions */
	RNG rng(seed);

	double* nodePr = new double[numNodes];
	Map<VectorXd> nodePrMap(nodePr, numNodes); /* use a map to access nodePr indirectly */
	nodePrMap.setConstant(1);

	BranchDistrib branch_dist;

	LocDistrib loc_dist(0, csLen - 1);
	if(!myLoci.empty()) /* restricted regions provided */
		loc_dist = LocDistrib(0, myLoci.size() - 1);

	SizeDistrib size_dist(meanSize, sdSize);

	GapDistrib gap_dist;

	double basePr[4] = {1, 1, 1, 1};
	Map<Vector4d> basePrMap(basePr, 4); /* use a map to access basePr indirectly */
	BaseDistrib base_dist(basePr);

	char rid[22]; // enough to hold all numbers up to 64-bits + a prefix char
	char desc[4096]; // enough to hold must descriptions
	const DegenAlphabet* abc = msa.getAbc();
	const PTUnrooted::ModelPtr& model = ptu.getModel();
	const Vector4d& pi = model->getPi();

	if(maxDist != EGriceLab::inf) { /* if -d specified */
		/* alter the node_dist weight */
		for(size_t i = 0; i < numNodes; ++i)
			if(ptu.getHeight(i) > maxDist)
				nodePrMap(i) = 0;
	}
	/* construct node dist w/ potentially modified weights */
	NodeDistrib node_dist(nodePr, nodePr + numNodes);

	/* generating random reads */
	if(mateFn.empty())
		infoLog << "Simulating single-end reads" << endl;
	else
		infoLog << "Simulating paired-end reads" << endl;
	for(long n = 1; n <= N;) {
		/* simulate a branch */
		size_t id = node_dist(rng);
		PTUnrooted::PTUNodePtr cNode = ptu.getNode(id);
		if(cNode->isRoot()) /* no parent branch available */
			continue;
		PTUnrooted::PTUNodePtr pNode = cNode->getParent();
		double v = ptu.getBranchLength(pNode, cNode);
		double rc = branch_dist(rng);
		if(ptu.getHeight(cNode) + v * rc > maxDist) /* this branching-point it too far from any leaf */
			continue;
		/* simulate a read range */
		int start, end, len;
		if(myLoci.empty()) { /* simulate from the entire CSLen */
			start = loc_dist(rng);
			len = size_dist(rng);
			if(len < minSize)
				len = minSize;
			if(maxSize > 0 && len > maxSize)
				len = maxSize;
			end = start + static_cast<int> (len);
			if(!(end < csLen)) /* outside consensus range */
				continue;
		}
		else { /* simulate from restricted regions */
			size_t i = loc_dist(rng);
			start = myLoci[i].csStart;
			end = myLoci[i].csEnd;
			len = end - start + 1;
		}

		/* simulate a read at [start, end] */
		sprintf(rid, "r%d", n);
//		PrimarySeq seq(abc, rid, "", desc);
		string seq;
		if(keepGap)
			seq.append(start, PAD_SYM);

		for(int j = start; j <= end; ++j) {
			bool isGap = gap_dist(rng) <= msa.gapWFrac(j);
			if(isGap) {
				if(keepGap)
					seq.push_back(GAP_SYM);
			}
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

		if(keepGap)
			seq.append(csLen - 1 - end, PAD_SYM);
		sprintf(desc, "ID=%ld->%ld;Name=\"%s\";AnnoDist=%f;csStart=%d;csEnd=%d;csLen=%d;InsertLen=%d;",
				cNode->getId(), pNode->getId(),
				rc < 0.5 ? cNode->getTaxon().c_str() : pNode->getTaxon().c_str(),
				rc < 0.5 ? v * rc : v * (1 - rc),
				start, end, end - start + 1, seq.length());

		/* output */
		PrimarySeq insert(abc, rid, seq, desc);
		seqOut.writeSeq(insert.trunc(0, readLen));
		if(mateOut.is_open())
			mateOut.writeSeq(insert.revcom().trunc(0, readLen));
		n++;
	}
}
