/*
 ============================================================================
 Name        : hmmufotu
 Author      : Qi Zheng
 Version     : v1.1
 Copyright   : Your copyright notice
 Description : Hello World in C++,
 ============================================================================
 */

#include <iostream>
#include <fstream>
#include <map>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <ctime>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include "HmmUFOtu.h"

using namespace std;
using namespace EGriceLab;
using namespace Eigen;

/** default values */
static const string ALPHABET = "dna";
static const double DEFAULT_MAX_PDIST = 0.03;
static const int DEFAULT_MIN_Q = 0;
static const int DEFAULT_SEED_LEN = 20;
static const int MAX_SEED_LEN = 25;
static const int MIN_SEED_LEN = 15;
static const int DEFAULT_SEED_REGION = 50;
static const double EST_PLACE_ERROR = 0.001;
static const string UNASSIGNED_TAXA = "Unknown";
static const double UNASSIGNED_LOGLIK = EGriceLab::nan;
static const string UNASSIGNED_ID = "NULL";
static const double UNASSIGNED_POSTQ = EGriceLab::nan;
static const double UNASSIGNED_DIST = EGriceLab::nan;
static const double UNASSIGNED_RATIO = EGriceLab::nan;

static const string PLACE_HEADER = "id\tdesc\tbranch_id\tbranch_ratio\tCS_start\tCS_end\ttaxa_anno\tanno_dist\tloglik\tpostQ";
static const string TAXA_SUMM_HEADER = "taxa_id\ttaxa_annotation\tcount";

typedef boost::unordered_set<PTUnrooted::PTUNodePtr> NodeSet;

/** struct to store placement information */
struct PTPlacement {
	/** construct a placement using given info */
	PTPlacement(const string& id, const string& taxa,
			double loglik, double ratio, double annoDist, double q = 0)
	: id(id), taxa(taxa), loglik(loglik), ratio(ratio), annoDist(annoDist), q(q)
	{ }

	static void calcQValues(vector<PTPlacement>& places);

	friend bool operator<(const PTPlacement& lhs, const PTPlacement& rhs);

	string id;
	string taxa;
	double loglik;
	double ratio;
	double annoDist;
	double q;
};

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Ultra-fast 16S read OTU assignment using profile-HMM and phylogenetic placement" << endl
		 << "Usage:    " << progName << "  <HmmUFOtu-DB> <READ-FILE1> [READ-FILE2] [-o OUTPUT] [options]" << endl
		 << "READ-FILE1  FILE               : sequence read file for the assembled/forward read" << endl
		 << "READ-FILE2  FILE               : sequence read file for the reverse read" << endl
		 << "Options:    -o  STR            : prefix for output files" << endl
		 << "            -L|--seed-len  INT : seed length used for banded-Hmm search [" << DEFAULT_SEED_LEN << "]" << endl
		 << "            -R  INT            : size of 5'/3' seed region for finding seeds [" << DEFAULT_SEED_REGION << "]" << endl
		 << "            -s  FLAG           : assume READ-FILE1 is single-end read instead of assembled read, if no READ-FILE2 provided" << endl
		 << "            -d|--pdist  DBL    : maximum p-dist between placed read and observed leaves during the optimal search [" << DEFAULT_MAX_PDIST << "]" << endl
		 << "            -q  DBL            : minimum Q-value (negative log10 of posterior placement probability) required for a read placement [" << DEFAULT_MIN_Q << "]" << endl
		 << "            -S|--seed  INT     : random seed used for banded HMM seed searches, for debug purpose" << endl
		 << "            -v  FLAG           : enable verbose information" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

vector<PTUnrooted::PTUNodePtr> getLeafHitsById(const set<unsigned>& idHits,
		const map<unsigned, PTUnrooted::PTUNodePtr>& id2leaf);

/** Align seq with hmm and csfm, returns the alignment or empty seq if failed */
PrimarySeq alignSeq(const BandedHMMP7& hmm, const CSFMIndex& csfm, const PrimarySeq& read,
		int seedLen, int seedRegion, BandedHMMP7::align_mode mode,
		int& csStart, int& csEnd);

/** Place a seq with known profile-HMM and Tree, modify the bestAnnotation accordingly */
PTPlacement placeSE(const PTUnrooted& ptu, const DigitalSeq& seq, int start, int end,
		double maxDist);

/** Place a mate of sequences with known profile-HMM and Tree, modify the bestAnnotation accordingly */
void placePE(const PTUnrooted& ptu, const DigitalSeq& fwdSeq, const DigitalSeq& revSeq,
		int fwdStart, int fwdEnd, int revStart, int revEnd,
		double maxDist, double minQ);

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, fwdFn, revFn, msaFn, csfmFn, hmmFn, ptuFn;
	string outName, outFn, alnFn;
	ifstream msaIn, csfmIn, hmmIn, ptuIn;
	string seqFmt;
	ofstream out;
	SeqIO alnOut;
	bool isAssembled = true; /* assume assembled seq if not paired-end */
	BandedHMMP7::align_mode mode;

	int seedLen = DEFAULT_SEED_LEN;
	int seedRegion = DEFAULT_SEED_REGION;
	double maxDist = DEFAULT_MAX_PDIST;
	double minQ = DEFAULT_MIN_Q;

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
		outName = cmdOpts.getOpt("-o");
	else {
		cerr << "Error: -o must be specified" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-L"))
		seedLen = ::atoi(cmdOpts.getOptStr("-L"));
	if(cmdOpts.hasOpt("--seed-len"))
		seedLen = ::atoi(cmdOpts.getOptStr("--seed-len"));
	if(!(MIN_SEED_LEN <= seedLen && seedLen <= MAX_SEED_LEN)) {
		cerr << "-L|--seed-len must be in range [" << MIN_SEED_LEN << ", " << MAX_SEED_LEN << "]" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-R"))
		seedRegion = ::atoi(cmdOpts.getOptStr("-R"));
	if(seedRegion < seedLen) {
		cerr << "-R cannot be smaller than -L" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-s"))
		isAssembled = false;

	if(cmdOpts.hasOpt("-d"))
		maxDist = ::atof(cmdOpts.getOptStr("-d"));
	if(cmdOpts.hasOpt("--pdist"))
		maxDist = ::atof(cmdOpts.getOptStr("--pdist"));

	if(cmdOpts.hasOpt("-q"))
		minQ = ::atof(cmdOpts.getOptStr("-q"));
	if(minQ < 0) {
		cerr << "-q must be non-negative" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-S"))
		seed = ::atoi(cmdOpts.getOptStr("-S"));
	if(cmdOpts.hasOpt("--seed"))
		seed = ::atoi(cmdOpts.getOptStr("--seed"));
	srand(seed);

	if(cmdOpts.hasOpt("-v"))
		ENABLE_INFO();

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

	/* set filenames */
	msaFn = dbName + MSA_FILE_SUFFIX;
	csfmFn = dbName + CSFM_FILE_SUFFIX;
	hmmFn = dbName + HMM_FILE_SUFFIX;
	ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;

	outFn = outName + "_placement.txt";
	alnFn = outName + "_aliged.fasta";

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

	/* open outputs */
	out.open(outFn.c_str());
	if(!out.is_open()) {
		cerr << "Unable to write to '" << outFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	alnOut.open(alnFn, "dna", "fasta", SeqIO::WRITE);
	if(!alnOut.is_open()) {
		cerr << "Unable to write to '" << alnFn << "': " << ::strerror(errno) << endl;
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

	CSFMIndex csfm;
	csfm.load(csfmIn);
	if(csfmIn.bad()) {
		cerr << "Failed to load CSFM-index '" << csfmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "CSFM-index loaded." << endl;
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
	infoLog << "HMM profile read." << endl;
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
	infoLog << "Phylogenetic tree loaded." << endl;

	const DegenAlphabet* abc = hmm.getNuclAbc();

	/* configure HMM mode */
	hmm.setSequenceMode(mode);
	hmm.wingRetract();

	/* get MSA index */
	const map<unsigned, PTUnrooted::PTUNodePtr>& id2leaf = ptu.getMSAIndex();

	/* process reads and output */
	out << PLACE_HEADER << endl;
	map<string, long> taxaCount;

	if(revFn.empty()) { /* single-ended */
		SeqIO seqIn(fwdFn, ALPHABET, seqFmt);
		while(seqIn.hasNext()) {
			const PrimarySeq& read = seqIn.nextSeq();
			cerr << "Aligning read: " << read.getId() << endl;
			int csStart = -1;
			int csEnd = -1;
			/* align the read */
			const PrimarySeq& aln = alignSeq(hmm, csfm, read, seedLen, seedRegion, mode, csStart, csEnd);
			if(aln.empty())
				continue;

			assert(aln.length() == csLen);
			DigitalSeq seq(aln);

			if(!alnFn.empty())
				alnOut.writeSeq(aln);

			/* place seq */
			const PTPlacement& place = placeSE(ptu, seq, csStart, csEnd, maxDist);

			if(::isnan(place.loglik))
				warningLog << "Unable to place read " << read.getId() << endl;

			taxaCount[place.taxa]++;

			/* write main output */
			out << read.getId() << "\t" << read.getDesc() << "\t"
				<< place.id << "\t" << place.ratio << "\t"
				<< csStart << "\t" << csEnd << "\t"
				<< place.taxa << "\t" << place.annoDist << "\t"
				<< place.loglik << "\t" << place.q << endl;
		}
	} /* end single-ended */
	else { /* paired-ended */

	} /* end paired-ended */

	return 0;
}

vector<PTUnrooted::PTUNodePtr> getLeafHitsById(const set<unsigned>& idHits,
		const map<unsigned, PTUnrooted::PTUNodePtr>& id2leaf) {
	vector<PTUnrooted::PTUNodePtr> leafHits;
	for(set<unsigned>::const_iterator id = idHits.begin(); id != idHits.end(); ++id)
		leafHits.push_back(id2leaf.at(*id));
	return leafHits;
}

PrimarySeq alignSeq(const BandedHMMP7& hmm, const CSFMIndex& csfm, const PrimarySeq& read,
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
	cerr << "nSeed: " << seqVpath.N << endl;

	/* banded HMM align */
	hmm.calcViterbiScores(seqVscore, seqVpath);
	float minCost = hmm.buildViterbiTrace(seqVscore, seqVpath);
	if(minCost == inf) {
		warningLog << "Warning: Unable to align read " << read.getId() << " , no Viterbi path found" << endl;
		return PrimarySeq(abc, read.getId(), "");
	}

	/* find seqStart and seqEnd */
	csStart = hmm.getCSLoc(seqVpath.alnStart) - 1;
	csEnd = hmm.getCSLoc(seqVpath.alnEnd) - 1;

	return hmm.buildGlobalAlignSeq(seqVscore, seqVpath);
}

PTPlacement placeSE(const PTUnrooted& ptu, const DigitalSeq& seq, int start, int end,
		double maxDist) {
	PTUnrooted::PTUNodePtr bestNode;
	/* Use a LAEP (leaf-ancestor-estimate-place) algorithm */
	const vector<PTUnrooted::PTUNodePtr>& leafHits = ptu.getLeafHitsByPDist(seq, maxDist, start, end);
	infoLog << "Found " << leafHits.size() << " leaf nodes for " << seq.getName() << endl;
	if(leafHits.empty())
		return PTPlacement(UNASSIGNED_ID, UNASSIGNED_TAXA, UNASSIGNED_LOGLIK, UNASSIGNED_RATIO, UNASSIGNED_DIST, UNASSIGNED_POSTQ);

	vector<PTUnrooted::PTUNodePtr> candidates;
	vector<PTPlacement> places;
	NodeSet nodeSeen;

	for(vector<PTUnrooted::PTUNodePtr>::const_iterator leaf = leafHits.begin(); leaf != leafHits.end(); ++leaf) {
		double prevLoglik = infV;
		for(PTUnrooted::PTUNodePtr node = *leaf; !node->isRoot() && nodeSeen.find(node) == nodeSeen.end(); node = node->getParent()) {
			nodeSeen.insert(node);
			/* estimate the placement */
			double loglik = ptu.estimateSeq(seq, node, node->getParent(), start, end);
			if(prevLoglik != infV && loglik - prevLoglik < prevLoglik * EST_PLACE_ERROR)
				break; /* stop search this lineage */
			candidates.push_back(node);
			prevLoglik = loglik;
		}
	}
	cerr << "Find " << candidates.size() << " branch candidates" << endl;

	for(vector<PTUnrooted::PTUNodePtr>::const_iterator node = candidates.begin(); node != candidates.end(); ++node) {
		/* place the read here */
		PTUnrooted subtree = ptu.copySubTree(*node, (*node)->getParent());

		const PTUnrooted::PTUNodePtr& v = subtree.getNode(0);
		const PTUnrooted::PTUNodePtr& u = subtree.getNode(1);
		double w0 = subtree.getBranchLength(u, v);

		double loglik = subtree.placeSeq(seq, u, v, start, end);
		const PTUnrooted::PTUNodePtr& r = subtree.getNode(2);
		const PTUnrooted::PTUNodePtr& n = subtree.getNode(3);

		/* calculate the node n annoDist */
		double wnr = subtree.getBranchLength(n, r);
		double wur = subtree.getBranchLength(u, r);
		double wvr = subtree.getBranchLength(v, r);
		double ratio = wur / w0;
		if(wvr <= wur) { /* r->v is shorter */
			n->setAnno(v->getAnno());
			n->setAnnoDist(v->getAnnoDist() + wvr + wnr);
		}
		else { /* u->r is shorter */
			n->setAnno(u->getAnno());
			if(u->getAnnoDist() == 0) /* self-annotated */
				n->setAnnoDist(wur + wnr);
			else /* annotated from v or ancestor */
				n->setAnnoDist(u->getAnnoDist() + (wnr - wur));
		}
		char branchID[64];
		sprintf(branchID, "%ld->%ld", (*node)->getId(), (*node)->getParent()->getId());
		places.push_back(PTPlacement(branchID, n->getTaxa(), loglik, ratio, n->getAnnoDist()));
	} /* end each candidate */

	cerr << "Tried " << places.size() << " places" << endl;
	/* sort placements descendingly */
	std::sort(places.rbegin(), places.rend());
	PTPlacement::calcQValues(places);

	return places.front();
}

void PTPlacement::calcQValues(vector<PTPlacement>& places) {
	if(places.empty())
		return;

	VectorXd loglik(places.size());
	VectorXd::Index i = 0;
	for(vector<PTPlacement>::const_iterator placement = places.begin(); placement != places.end(); ++placement)
		loglik(i++) = placement->loglik;

	/* scale and normalize */
	VectorXd p = (loglik.array() - loglik.maxCoeff()).exp();
	p /= p.sum();

	/* assign q-values */
	for(vector<PTPlacement>::size_type i = 0; i < places.size(); ++i)
		places[i].q = -::log10(1 - p(i));
}

inline bool operator<(const PTPlacement& lhs, const PTPlacement& rhs) {
	return lhs.loglik < rhs.loglik;
}
