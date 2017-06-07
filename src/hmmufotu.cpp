/*
 ============================================================================
 Name        : hmmufotu
 Author      : Qi Zheng
 Version     : v1.1
 Copyright   : Your copyright notice
 Description : Main program of the HmmUFOtu project
 ============================================================================
 */

#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <map>
#include <algorithm>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "HmmUFOtu.h"

using namespace std;
using namespace EGriceLab;
using namespace Eigen;

/* default values */
static const double DEFAULT_MAX_PDIST = inf;
static const size_t DEFAULT_MAX_LOCATION = 50;
static const int DEFAULT_MIN_Q = 0;
static const int MAX_Q = 250; /* maximum allowed Q value */
static const int DEFAULT_SEED_LEN = 20;
static const int MAX_SEED_LEN = 25;
static const int MIN_SEED_LEN = 15;
static const int DEFAULT_SEED_REGION = 50;
static const double DEFAULT_MAX_PLACE_ERROR = 20;
static const int DEFAULT_NUM_THREADS = 1;
static const long UNASSIGNED_TAXONID = -1;
static const string UNASSIGNED_TAXONNAME = "UNASSIGNED";
static const double UNASSIGNED_LOGLIK = EGriceLab::nan;
static const string UNASSIGNED_ID = "NULL";
static const double UNASSIGNED_POSTQ = EGriceLab::nan;
static const double UNASSIGNED_DIST = EGriceLab::nan;
static const double UNASSIGNED_RATIO = EGriceLab::nan;

static const string PLACE_HEADER = "id\tdescription\tCS_start\tCS_end\talignment\tbranch_id\tbranch_ratio\ttaxon_id\ttaxon_anno\tanno_dist\tloglik\tQ_placement\tQ_taxon";

/** prior probability method */
enum PRIOR_TYPE {
	UNIFORM,
	HEIGHT
};

/** struct to store potential placement location information reads */
struct PTLoc {
	/** construct a location using given info */
	PTLoc(const PTUnrooted::PTUNodePtr& node, double dist)
	: node(node), dist(dist)
	{ }

	friend bool operator<(const PTLoc& lhs, const PTLoc& rhs);

	PTUnrooted::PTUNodePtr node;
	double dist;
};

/** struct to store placement information */
struct PTPlacement {
//	/** default constructor */
	PTPlacement() { }

	/** construct a placement with basic info and optionally auxilary info */
	PTPlacement(const PTUnrooted::PTUNodePtr& cNode, const PTUnrooted::PTUNodePtr& pNode,
			double ratio, double wnr, double loglik,
			double height = 0, double annoDist = 0,
			double qPlace = 0, double qTaxonomy = 0)
	: cNode(cNode), pNode(pNode), ratio(ratio), wnr(wnr), loglik(loglik),
	  height(height), annoDist(annoDist), qPlace(qPlace), qTaxon(qTaxonomy)
	{ }

	long getTaxonId() const {
		return ratio <= 0.5 ? cNode->getId() : pNode->getId();
	}

	string getTaxonName() const {
		return ratio <= 0.5 ? cNode->getAnno() : pNode->getAnno();
	}

	string getId() const {
		return boost::lexical_cast<string> (cNode->getId()) + "->" + boost::lexical_cast<string> (cNode->getParent()->getId());
	}

	static void calcQValues(vector<PTPlacement>& places, PRIOR_TYPE type);

	friend bool compareByLoglik(const PTPlacement& lhs, const PTPlacement& rhs);
	friend bool compareByQTaxon(const PTPlacement& lhs, const PTPlacement& rhs);
	friend bool compareByQPlace(const PTPlacement& lhs, const PTPlacement& rhs);
	friend ostream& operator<<(ostream& out, const PTPlacement& place);

	PTUnrooted::PTUNodePtr cNode;
	PTUnrooted::PTUNodePtr pNode;
	double ratio; /* placement ratio */
	double wnr;   /* new branch length */
	double loglik;
	double annoDist;
	double height;
	double qPlace;
	double qTaxon;
};

inline bool compareByLoglik(const PTPlacement& lhs, const PTPlacement& rhs) {
	return lhs.loglik < rhs.loglik;
}

inline bool compareByQPlace(const PTPlacement& lhs, const PTPlacement& rhs) {
	return lhs.qPlace < rhs.qPlace;
}

inline bool compareByQTaxon(const PTPlacement& lhs, const PTPlacement& rhs) {
	return lhs.qTaxon < rhs.qTaxon;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Ultra-fast 16S read OTU assignment using profile-HMM and phylogenetic placement" << endl
		 << "Usage:    " << progName << "  <HmmUFOtu-DB> <READ-FILE1> [READ-FILE2] [options]" << endl
		 << "READ-FILE1  FILE               : sequence read file for the assembled/forward read" << endl
		 << "READ-FILE2  FILE               : sequence read file for the reverse read" << endl
		 << "Options:    -o  FILE           : write the PLACEMENT output to FILE instead of stdout" << endl
		 << "            -a  FILE           : in addition to the placement output, write the read alignment to FILE" << endl
		 << "            -L|--seed-len  INT : seed length used for banded-Hmm search [" << DEFAULT_SEED_LEN << "]" << endl
		 << "            -R  INT            : size of 5'/3' seed region for finding seeds [" << DEFAULT_SEED_REGION << "]" << endl
		 << "            -s  FLAG           : assume READ-FILE1 is single-end read instead of assembled read, if no READ-FILE2 provided" << endl
		 << "            -d|--pdist  DBL    : maximum p-dist between read and tree nodes during the seed search stage [" << DEFAULT_MAX_PDIST << "]" << endl
		 << "            -N  INT            : max number of most potential locations (based on p-dist) to try initial estimation [" << DEFAULT_MAX_LOCATION << "]" << endl
		 << "            -e  DBL            : max placement error between fast estimation and accurate placement [" << DEFAULT_MAX_PLACE_ERROR << "]" << endl
		 << "            --ML  FLAG         : use maximum likelihood for placement, do not calculate posterior p-values, this will ignore -q and --prior options" << endl
		 << "            -q  INT            : minimum Q-value (negative log10 of posterior placement probability) required for a read placement [" << DEFAULT_MIN_Q << "]" << endl
		 << "            --prior  STR       : method for calculating prior probability of a placement, either 'uniform' (uniform prior) or 'height' (rooted distance to leaves)" << endl
		 << "            -S|--seed  INT     : random seed used for banded HMM seed searches, for debug purpose" << endl
#ifdef _OPENMP
		 << "            -p|--process INT   : number of threads/cpus used for parallel processing" << endl
#endif
		 << "            -v  FLAG           : enable verbose information, you may set multiple -v for more details" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

/** Align seq with hmm and csfm, returns the alignment and update csStart and csEnd */
string alignSeq(const BandedHMMP7& hmm, const CSFMIndex& csfm, const PrimarySeq& read,
		int seedLen, int seedRegion, BandedHMMP7::align_mode mode,
		int& csStart, int& csEnd);

/**
 * Get seed placement locations by checking p-dist between a given seq and observed/inferred seq of nodes
 * @param ptu  PTUnrooted
 * @param seq  sequence to be placed
 * @param start  0-based start
 * @param end  0-based end
 * @param maxDist  maximum p-Distance
 * */
vector<PTLoc> getSeed(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, double maxDist);

/** Get estimated placement for a seq at given locations */
vector<PTPlacement> estimateSeq(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, const vector<PTLoc>& locs);

/** Get accurate placement for a seq given the estimated placements */
vector<PTPlacement>& placeSeq(const PTUnrooted& ptu, const DigitalSeq& seq, int start, int end,
		vector<PTPlacement>& places);

/**
 * calculate prior probability at log-scale
 * @param place  a placement
 * @param type  prior type
 * @param h  base height of this placement (for cNode)
 * @return  log prior always no greater than 0
 */
double logPriorPr(const PTPlacement& place, PRIOR_TYPE type, double h) {
	double logP;
	switch(type) {
	case UNIFORM:
		logP = -0;
		break;
	case HEIGHT:
		logP = -(place.annoDist - place.wnr + h);
		break;
	}
	return logP;
}

/** calculate prior proability of a placement using given type of method */
double priorPr(const PTPlacement& place, PRIOR_TYPE type, double h) {
	return ::exp(logPriorPr(place, type, h));
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, fwdFn, revFn, msaFn, csfmFn, hmmFn, ptuFn;
	string outFn, alnFn;
	ifstream msaIn, csfmIn, hmmIn, ptuIn;
	string seqFmt;
	ofstream of;
	SeqIO fwdIn, revIn, alnOut;
	bool isAssembled = true; /* assume assembled seq if not paired-end */
	BandedHMMP7::align_mode mode;

	int seedLen = DEFAULT_SEED_LEN;
	int seedRegion = DEFAULT_SEED_REGION;
	double maxDist = DEFAULT_MAX_PDIST;
	int maxLocs = DEFAULT_MAX_LOCATION;
	double maxError = DEFAULT_MAX_PLACE_ERROR;
	bool onlyML = false;
	PRIOR_TYPE myPrior = UNIFORM;
	int minQ = DEFAULT_MIN_Q;
	int nThreads = DEFAULT_NUM_THREADS;

	unsigned seed = time(NULL); // using time as default seed

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(!(cmdOpts.numMainOpts() == 2 || cmdOpts.numMainOpts() == 3)) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	dbName = cmdOpts.getMainOpt(0);
	fwdFn = cmdOpts.getMainOpt(1);
	if(cmdOpts.numMainOpts() == 3)
		revFn = cmdOpts.getMainOpt(2);

	if(cmdOpts.hasOpt("-o"))
		outFn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("-a"))
		alnFn = cmdOpts.getOpt("-a");

	if(cmdOpts.hasOpt("-L"))
		seedLen = ::atoi(cmdOpts.getOptStr("-L"));
	if(cmdOpts.hasOpt("--seed-len"))
		seedLen = ::atoi(cmdOpts.getOptStr("--seed-len"));

	if(cmdOpts.hasOpt("-R"))
		seedRegion = ::atoi(cmdOpts.getOptStr("-R"));

	if(cmdOpts.hasOpt("-s"))
		isAssembled = false;

	if(cmdOpts.hasOpt("-d"))
		maxDist = ::atof(cmdOpts.getOptStr("-d"));
	if(cmdOpts.hasOpt("--pdist"))
		maxDist = ::atof(cmdOpts.getOptStr("--pdist"));

	if(cmdOpts.hasOpt("-N"))
		maxLocs = ::atoi(cmdOpts.getOptStr("-N"));

	if(cmdOpts.hasOpt("-e"))
		maxError = ::atof(cmdOpts.getOptStr("-e"));

	if(cmdOpts.hasOpt("--ML"))
		onlyML = true;

	if(cmdOpts.hasOpt("--prior")) {
		if(cmdOpts.getOpt("--prior") == "uniform")
			myPrior = UNIFORM;
		else if(cmdOpts.getOpt("--prior") == "height")
			myPrior = HEIGHT;
		else {
			cerr << "Unsupported prior specified, check the --prior option" << endl;
			return EXIT_FAILURE;
		}
	}

	if(cmdOpts.hasOpt("-q"))
		minQ = ::atoi(cmdOpts.getOptStr("-q"));

	if(cmdOpts.hasOpt("-S"))
		seed = ::atoi(cmdOpts.getOptStr("-S"));
	if(cmdOpts.hasOpt("--seed"))
		seed = ::atoi(cmdOpts.getOptStr("--seed"));
	srand(seed);

#ifdef _OPENMP
	if(cmdOpts.hasOpt("-p"))
		nThreads = ::atoi(cmdOpts.getOptStr("-p"));
	if(cmdOpts.hasOpt("--process"))
		nThreads = ::atoi(cmdOpts.getOptStr("--process"));
#endif

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* guess seq format */
	if(StringUtils::endsWith(fwdFn, ".fasta") || StringUtils::endsWith(fwdFn, ".fas")
		|| StringUtils::endsWith(fwdFn, ".fa") || StringUtils::endsWith(fwdFn, ".fna"))
		seqFmt = "fasta";
	else if(StringUtils::endsWith(fwdFn, ".fastq") || StringUtils::endsWith(fwdFn, ".fq"))
		seqFmt = "fastq";
	else {
		cerr << "Unrecognized format of MSA file '" << fwdFn << "'" << endl;
		return EXIT_FAILURE;
	}

	/* validate options */
	if(!(MIN_SEED_LEN <= seedLen && seedLen <= MAX_SEED_LEN)) {
		cerr << "-L|--seed-len must be in range [" << MIN_SEED_LEN << ", " << MAX_SEED_LEN << "]" << endl;
		return EXIT_FAILURE;
	}
	if(seedRegion < seedLen) {
		cerr << "-R cannot be smaller than -L" << endl;
		return EXIT_FAILURE;
	}
	if(minQ < 0) {
		cerr << "-q must be non-negative" << endl;
		return EXIT_FAILURE;
	}
	if(!(maxDist > 0)) {
		cerr << "-d must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(!(maxLocs > 0)) {
		cerr << "-N must be positive" << endl;
		return EXIT_FAILURE;
	}
	if(!(maxError > 0)) {
		cerr << "-e must be positive" << endl;
		return EXIT_FAILURE;
	}
#ifdef _OPENMP
	if(!(nThreads > 0)) {
		cerr << "-p|--process must be positive" << endl;
		return EXIT_FAILURE;
	}
	omp_set_num_threads(nThreads);
#endif

	if(onlyML)
		minQ = 0;

	/* set filenames */
	msaFn = dbName + MSA_FILE_SUFFIX;
	csfmFn = dbName + CSFM_FILE_SUFFIX;
	hmmFn = dbName + HMM_FILE_SUFFIX;
	ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;

	/* set HMM align mode */
	mode = !revFn.empty() /* paired-end */ || isAssembled ? BandedHMMP7::GLOBAL : BandedHMMP7::NGCL;

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

	fwdIn.open(fwdFn, AlphabetFactory::nuclAbc, seqFmt);
	if(!fwdIn.is_open()) {
		cerr << "Unable to open seq file '" << fwdFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(!revFn.empty()) {
		revIn.open(revFn, AlphabetFactory::nuclAbc, seqFmt);
		if(!revIn.is_open()) {
			cerr << "Unable to open mate file '" << revFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* open outputs */
	if(!outFn.empty()) {
		of.open(outFn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to '" << outFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}
	ostream& out = of.is_open() ? of : cout;

	if(!alnFn.empty()) {
		alnOut.open(alnFn, "dna", "fasta", SeqIO::WRITE);
		if(!alnOut.is_open()) {
			cerr << "Unable to write to '" << alnFn << "': " << ::strerror(errno) << endl;
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

	CSFMIndex csfm;
	csfm.load(csfmIn);
	if(csfmIn.bad()) {
		cerr << "Failed to load CSFM-index '" << csfmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "CSFM-index loaded" << endl;
	if(csfm.getCSLen() != csLen) {
		cerr << "Error: Unmatched CS length between CSFM-index and MSA data" << endl;
		return EXIT_FAILURE;
	}

	BandedHMMP7 hmm;
	hmmIn >> hmm;
	if(hmmIn.bad()) {
		cerr << "Unable to read HMM profile '" << hmmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "HMM profile read" << endl;
	if(hmm.getProfileSize() > csLen) {
		cerr << "Error: HMM profile size is found greater than the MSA CS length" << endl;
		return EXIT_FAILURE;
	}

	PTUnrooted ptu;
	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Unable to load Phylogenetic tree data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Phylogenetic tree loaded" << endl;

	const DegenAlphabet* abc = hmm.getNuclAbc();

	/* configure HMM mode */
	hmm.setSequenceMode(mode);
	hmm.wingRetract();


	/* process reads and output */
	out << PLACE_HEADER << endl;
#pragma omp parallel shared(csfm, hmm, ptu, fwdIn, revIn)
#pragma omp single
	{
		while(fwdIn.hasNext() && (!revIn.is_open() || revIn.hasNext())) {
			string id;
			string desc;
			PrimarySeq fwdRead, revRead;
			fwdRead = fwdIn.nextSeq();
			id = fwdRead.getId();
			desc = fwdRead.getDesc();
			if(revIn.is_open()) { /* paired-ended */
				revRead = revIn.nextSeq().revcom();
				if(fwdRead.getId() != revRead.getId())
					warningLog << "Warning: Paired-end read doesn't match at " << fwdRead.getId() << " vs. " << revRead.getId() << endl;
			}

#pragma omp task
			{
				int csStart = 0; /* 1-based consensus start */
				int csEnd = 0;   /* 1-based consensus end */
				string aln;
				PTPlacement bestPlace;

				/* align fwdRead */
				aln = alignSeq(hmm, csfm, fwdRead, seedLen, seedRegion, mode, csStart, csEnd);
				assert(aln.length() == csLen);
				if(revIn.is_open()) { /* align revRead */
					// cerr << "Aligning mate: " << fwdRead.getId() << endl;
					int revStart = 0;
					int revEnd = 0;
					string revAln = alignSeq(hmm, csfm, revRead, seedLen, seedRegion, mode, revStart, revEnd);
					assert(revAln.length() == csLen);
					if(!(csStart <= revStart && csEnd <= revEnd))
#pragma omp critical(warning)
						warningLog << "Warning: Potential incorrectly oriented read found at " << id << endl;

					/* update alignment */
					csStart = ::min(csStart, revStart);
					csEnd = ::max(csEnd, revEnd);
					BandedHMMP7::mergeWith(aln, revAln);
				}
				assert(1 <= csStart && csStart <= csEnd && csEnd <= csLen);

				if(alnOut.is_open()) { /* write the alignment seq to output */
					string desc = fwdRead.getDesc();
					desc += ";csStart=" + boost::lexical_cast<string>(csStart) + ";csEnd=" + boost::lexical_cast<string>(csEnd) + ";";
#pragma omp critical(writeAln)
					alnOut.writeSeq(PrimarySeq(abc, id, aln, fwdRead.getDesc()));
				}

				DigitalSeq seq(abc, id, aln);
				/* place seq with seed-estimate-place (SEP) algorithm */
				/* get potential locs */
				vector<PTLoc> locs = getSeed(ptu, seq, csStart - 1, csEnd - 1, maxDist);
				if(locs.empty()) {
#pragma omp critical(warning)
					warningLog << "Warning: No seed loci found at " << id << endl;
				}
				else {
					std::sort(locs.begin(), locs.end()); /* sort by dist */
					if(locs.size() > maxLocs)
						locs.erase(locs.end() - (locs.size() - maxLocs), locs.end()); /* remove last maxLocs elements */
					//	cerr << "Found " << locs.size() << " potential placement locations" << endl;

					/* estimate placement */
					vector<PTPlacement> places = estimateSeq(ptu, seq, csStart - 1, csEnd - 1, locs);
					std::sort(places.rbegin(), places.rend(), ::compareByLoglik); /* sort places decently by estimated loglik */
					double bestEstLoglik = places.front().loglik;
					vector<PTPlacement>::const_iterator it;
					for(it = places.begin(); it != places.end(); ++it) {
						if(::abs(it->loglik - bestEstLoglik) > maxError)
							break;
					}
					places.erase(it, places.end()); /* remove too bad placements */

					/* accurate placement */
					placeSeq(ptu, seq, csStart - 1, csEnd - 1, places);
					if(onlyML) { /* don't calculate q-values */
						std::sort(places.rbegin(), places.rend(), ::compareByLoglik); /* sort places decently by real loglik */
					}
					else { /* calculate q-values */
						PTPlacement::calcQValues(places, myPrior);
						std::sort(places.rbegin(), places.rend(), ::compareByQTaxon); /* sort places decently by QTaxon */
					}

					bestPlace = places.front();

				/* write main output */
#pragma omp critical(writePlace)
				out << id << "\t" << desc << "\t" << csStart << "\t" << csEnd << "\t" << aln << "\t"
						<< bestPlace << endl;
				}
			} /* end task */
		} /* end each read/pair */
#pragma omp taskwait
	} /* end single */
	return 0;
}

string alignSeq(const BandedHMMP7& hmm, const CSFMIndex& csfm, const PrimarySeq& read,
		int seedLen, int seedRegion, BandedHMMP7::align_mode mode,
		int& csStart, int& csEnd) {
	const DegenAlphabet* abc = hmm.getNuclAbc();

	BandedHMMP7::ViterbiScores seqVscore = hmm.initViterbiScores(read); // construct an empty reusable score
	BandedHMMP7::ViterbiAlignPath seqVpath = hmm.initViterbiAlignPath(read.length()); // construct an empty reusable path

	int regionLen = seedRegion < read.length() ? seedRegion : read.length(); /* search region */
	/* find seed in 5' */
	for(int seedFrom = 0; seedFrom < regionLen - seedLen + 1; ++seedFrom) {
		PrimarySeq seed(abc, read.getId(), read.subseq(seedFrom, seedLen));
		const CSLoc& loc = csfm.locateOne(seed.getSeq());
		if(loc.start > 0 && loc.end > 0) /* a read seed located */ {
//					fprintf(stderr, "start:%d end:%d from:%d to:%d  CSLen:%d CS:%s\n", loc.start, loc.end, seedFrom + 1, seedFrom + seedLen, loc.CS.length(), loc.CS.c_str());
			hmm.addKnownAlignPath(seqVpath, loc, seedFrom + 1, seedFrom + seedLen); /* seed_from and seed_to are 1-based */
			const set<unsigned>& hits = csfm.locateIndex(seed.getSeq());
			break; /* only one 5'-seed necessary */
		}
	}
	/* find seed in 3', if requested */
	if(mode == BandedHMMP7::GLOBAL) {
		for(int seedTo = read.length() - 1; seedTo - seedLen + 1 >= read.length() - regionLen; --seedTo) {
			PrimarySeq seed(abc, read.getId(), read.subseq(seedTo - seedLen + 1, seedLen));
			const CSLoc& loc = csfm.locateOne(seed.getSeq());
			if(loc.start > 0 && loc.end > 0) /* a read seed located */ {
//						fprintf(stderr, "start:%d end:%d from:%d to:%d  CSLen:%d CS:%s\n", loc.start, loc.end, seedTo - seedLen + 2, seedTo + 1, loc.CS.length(), loc.CS.c_str());
				hmm.addKnownAlignPath(seqVpath, loc, seedTo - seedLen + 2, seedTo + 1); /* seed_from and seed_to are 1-based */
				const set<unsigned>& hits = csfm.locateIndex(seed.getSeq());
				break; /* only one 5'-seed necessary */
			}
		}
	}
//	cerr << "nSeed: " << seqVpath.N << endl;

	/* banded HMM align */
	if(seqVpath.N > 0) { /* use banded Viterbi algorithm */
		hmm.calcViterbiScores(seqVscore, seqVpath);
		if(seqVscore.S.minCoeff() == inf) { /* banded version failed */
			debugLog << "Banded HMM algorithm didn't find a potential Viterbi path, returning to regular HMM" << endl;
			hmm.resetViterbiScores(seqVscore, read);
//			hmm.resetViterbiAlignPath(seqVpath, L);
			hmm.calcViterbiScores(seqVscore);
		}
	}
	else
		hmm.calcViterbiScores(seqVscore); /* use original Viterbi algorithm */
	float minCost = hmm.buildViterbiTrace(seqVscore, seqVpath);

	assert(minCost != inf);

	/* find seqStart and seqEnd */
	csStart = hmm.getCSLoc(seqVpath.alnStart);
	csEnd = hmm.getCSLoc(seqVpath.alnEnd);

	return hmm.buildGlobalAlign(seqVscore, seqVpath);
}

vector<PTLoc> getSeed(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, double maxDist) {
	vector<PTLoc> locs; /* candidate locations */
	/* get potential placement locations based on pDist to observed or inferred sequences */
	for(vector<PTUnrooted::PTUNodePtr>::size_type i = 0; i < ptu.numNodes(); ++i) {
		PTUnrooted::PTUNodePtr node = ptu.getNode(i);
		if(node->isRoot())
			continue;
		double pDist = SeqUtils::pDist(node->getSeq(), seq, start, end);
		if(pDist <= maxDist)
			locs.push_back(PTLoc(node, pDist));
	}
	return locs;
}

vector<PTPlacement> estimateSeq(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, const vector<PTLoc>& locs) {
	vector<PTPlacement> places;
	for(vector<PTLoc>::const_iterator loc = locs.begin(); loc < locs.end(); ++loc) {
		const PTUnrooted::PTUNodePtr& cNode = loc->node;
		const PTUnrooted::PTUNodePtr& pNode = cNode->getParent();
		double cDist = loc->dist;
		double pDist = SeqUtils::pDist(pNode->getSeq(), seq, start, end);
		double ratio = cDist / (cDist + pDist);
		if(::isnan(ratio)) // unable to estimate the ratio
			ratio = 0.5;
//		cerr << "Estimating at " << node->getId() << " cDist: " << cDist << " pDist: " << pDist << " ratio: " << ratio << endl;
		/* estimate the placement */
		double wnr;
		double loglik = ptu.estimateSeq(seq, cNode, pNode, start, end, ratio, wnr);
		places.push_back(PTPlacement(cNode, pNode, ratio, wnr, loglik));
	}
	return places;
}

/** Get accurate placement for a seq given the estimated placements */
vector<PTPlacement>& placeSeq(const PTUnrooted& ptu, const DigitalSeq& seq, int start, int end,
		vector<PTPlacement>& places) {
	/* accurate placement using estimated values */
	for(vector<PTPlacement>::iterator place = places.begin(); place != places.end(); ++place) {
		double ratio0 = place->ratio;
		double wnr0 = place->wnr;
		double loglik0 = place->loglik;

//		cerr << "Estimated placement ratio0: " << ratio0 << " wnr: " << wnr0 << " loglik: " << loglik0 << endl;
		PTUnrooted subtree = ptu.copySubTree(place->cNode, place->pNode);
		const PTUnrooted::PTUNodePtr& v = subtree.getNode(0);
		const PTUnrooted::PTUNodePtr& u = subtree.getNode(1);
		double w0 = subtree.getBranchLength(u, v);

		double loglik = subtree.placeSeq(seq, u, v, start, end, ratio0, wnr0);
		const PTUnrooted::PTUNodePtr& r = subtree.getNode(2);
		const PTUnrooted::PTUNodePtr& n = subtree.getNode(3);

		/* update placement info */
		double wnr = subtree.getBranchLength(n, r);
		double wur = subtree.getBranchLength(u, r);
		double wvr = w0 - wur;
		place->ratio = wur / w0;
//		cerr << "delta loglik: " << (loglik - loglik0) << endl;
		place->height = ptu.getHeight(place->cNode) + wur;
		place->annoDist = wvr <= wur ? wvr + wnr : wur + wnr;
		/* update other placement info */
		place->wnr = wnr;
		place->loglik = loglik;
	} /* end each candidate placement */

	return places;
}

void PTPlacement::calcQValues(vector<PTPlacement>& places, PRIOR_TYPE type) {
	if(places.empty())
		return;

	/* explore all placements */
	VectorXd ppPlace(places.size()); /* posterior logP at placement */
	map<string, double> ppTaxon; /* posterior logP at taxon */
	double ppTaxNorm = infV; /* log(0) */

	VectorXd::Index i = 0;
	for(vector<PTPlacement>::const_iterator placement = places.begin(); placement != places.end(); ++placement) {
		double p = placement->loglik + logPriorPr(*placement, type, placement->height);
		ppPlace(i++) = p;
		string taxonomy = placement->getTaxonName();
		if(ppTaxon.find(taxonomy) == ppTaxon.end())
			ppTaxon[taxonomy] = p;
		else
			ppTaxon[taxonomy] = EGriceLab::Math::add_scaled(ppTaxon[taxonomy], p);
		ppTaxNorm = EGriceLab::Math::add_scaled(ppTaxNorm, p);
	}
	/* scale and normalize llPlace */
	VectorXd p = (ppPlace.array() - ppPlace.maxCoeff()).exp();
	p /= p.sum();
	/* calculate qPlace */
	for(vector<PTPlacement>::size_type i = 0; i < places.size(); ++i) {
		double q = EGriceLab::Math::p2q(1 - p(i));
		places[i].qPlace = q > MAX_Q ? MAX_Q : q;
	}

	/* calculate qTaxonomy */
	for(vector<PTPlacement>::iterator placement = places.begin(); placement != places.end(); ++placement) {
		double q = EGriceLab::Math::p2q(1 - ::exp(ppTaxon[placement->getTaxonName()] - ppTaxNorm));
		placement->qTaxon = q > MAX_Q ? MAX_Q : q;
	}
}

inline bool operator<(const PTLoc& lhs, const PTLoc& rhs) {
	return lhs.dist < rhs.dist;
}


inline ostream& operator<<(ostream& out, const PTPlacement& place) {
	if(place.cNode != NULL && place.pNode != NULL)
		out << place.getId() << "\t" << place.ratio << "\t"
		<< place.getTaxonId() << "\t" << place.getTaxonName() << "\t"
		<< place.annoDist << "\t" << place.loglik << "\t" << place.qPlace << "\t" << place.qTaxon;
	else
		out << UNASSIGNED_ID << "\t" << UNASSIGNED_RATIO << "\t"
		<< UNASSIGNED_TAXONID << "\t" << UNASSIGNED_TAXONNAME << "\t"
		<< UNASSIGNED_DIST << "\t" << UNASSIGNED_LOGLIK << "\t" << UNASSIGNED_POSTQ << "\t" << UNASSIGNED_POSTQ;
	return out;
}

