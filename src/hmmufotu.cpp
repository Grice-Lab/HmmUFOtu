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

/** default values */
static const string ALPHABET = "dna";
static const double DEFAULT_MAX_PDIST = 0.03;
static const int DEFAULT_MIN_Q = 0;
static const int DEFAULT_SEED_LEN = 20;
static const int MAX_SEED_LEN = 25;
static const int MIN_SEED_LEN = 15;
static const int DEFAULT_SEED_REGION = 50;
static const string PLACE_HEADER = "id\tdesc\tbranch_id\tbranch_name\tplace_region\tphylogenetic_annotation";
static const string SUMM_HEADER = "taxa\tcount";

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Ultra-fast 16S read OTU assignment using profile-HMM and phylogenetic placement" << endl
		 << "Usage:    " << progName << "  <HmmUFOtu-DB> <READ-FILE1> [READ-FILE2] [options]" << endl
		 << "READ-FILE1  FILE               : sequence read file for the assembled/forward read" << endl
		 << "READ-FILE2  FILE               : sequence read file for the reverse read" << endl
		 << "Options:    -o  FILE           : write output to FILE instead of stdout" << endl
		 << "            -L|--seed-len  INT : seed length used for banded-Hmm search [" << DEFAULT_SEED_LEN << "]" << endl
		 << "            -R  INT            : size of 5'/3' seed region for finding seeds [" << DEFAULT_SEED_REGION << "]" << endl
		 << "            -s  FLAG           : assume READ-FILE1 is single-end read instead of assembled read, if no READ-FILE2 provided" << endl
		 << "            -T  FILE           : in addition to the placement output, write the taxa summary table to TAXA-FILE" << endl
		 << "            -d|--pdist  DBL    : maximum p-dist between placed read and observed leaves during the optimal search [" << DEFAULT_MAX_PDIST << "]" << endl
		 << "            -q  DBL            : minimum Q-value (negative log10 of the error) required for a read placement [" << DEFAULT_MIN_Q << "]" << endl
		 << "            -S|--seed  INT     : random seed used for banded HMM seed searches, for debug purpose" << endl
		 << "            -v  FLAG           : enable verbose information" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

vector<PTUnrooted::PTUNodePtr> getLeafHitsById(const set<unsigned>& idHits,
		const map<unsigned, PTUnrooted::PTUNodePtr>& id2leaf);

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, fwdFn, revFn, msaFn, csfmFn, hmmFn, ptuFn;
	string outFn, taxaFn;
	ifstream msaIn, csfmIn, hmmIn, ptuIn;
	string seqFmt;
	ofstream of, taxaOut;
	bool isAssembled = true; /* assume assembled seq if not paired-end */
	BandedHMMP7::align_mode mode;

	int seedLen = DEFAULT_SEED_LEN;
	int seedRegion = DEFAULT_SEED_REGION;
	double maxDist = DEFAULT_MAX_PDIST;
	double minQ = DEFAULT_MIN_Q;

	unsigned seed = time(NULL); // using time as default seed

	typedef boost::unordered_set<PTUnrooted::PTUNodePtr> NodeSet;

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
	if(!outFn.empty()) {
		of.open(outFn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to '" << outFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}
	ostream& out = of.is_open() ? of : cout;

	if(!taxaFn.empty()) {
		taxaOut.open(taxaFn.c_str());
		if(!taxaOut.is_open()) {
			cerr << "Unable to write to '" << taxaFn << "': " << ::strerror(errno) << endl;
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

	/* process reads */
	if(revFn.empty()) { /* single-ended */
		SeqIO seqIn(fwdFn, ALPHABET, seqFmt);
		while(seqIn.hasNext()) {
			const PrimarySeq& read = seqIn.nextSeq();
			if(mode == BandedHMMP7::GLOBAL && read.length() > csLen) {
				warningLog << "Warning: read length cannot longer than consensus length for assmbled read mode for read "
						<< read.getId() << " length " << read.length() << " CSLen: " << csLen << endl;
				continue;
			}
			BandedHMMP7::ViterbiScores seqVscore = hmm.initViterbiScores(read); // construct an empty reusable score
			BandedHMMP7::ViterbiAlignPath seqVpath = hmm.initViterbiAlignPath(read.length()); // construct an empty reusable path
			set<unsigned> seqIdHits;

			int regionLen = seedRegion < read.length() ? seedRegion : read.length(); /* search region */
			/* find seed in 5' */
			for(int seedFrom = 0; seedFrom < regionLen - seedLen + 1; ++seedFrom) {
				PrimarySeq seed(abc, read.getId(), read.subseq(seedFrom, seedLen));
				const CSLoc& loc = csfm.locateOne(seed.getSeq());
				if(loc.start > 0 && loc.end > 0) /* a read seed located */ {
//					fprintf(stderr, "start:%d end:%d from:%d to:%d  CSLen:%d CS:%s\n", loc.start, loc.end, seedFrom + 1, seedFrom + seedLen, loc.CS.length(), loc.CS.c_str());
					hmm.addKnownAlignPath(seqVpath, loc, seedFrom + 1, seedFrom + seedLen); /* seed_from and seed_to are 1-based */
					const set<unsigned>& hits = csfm.locateIndex(seed.getSeq());
					seqIdHits.insert(hits.begin(), hits.end());
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
						seqIdHits.insert(hits.begin(), hits.end());
						break; /* only one 5'-seed necessary */
					}
				}
			}
			cerr << "nSeed: " << seqVpath.N << " idxHits: " << seqIdHits.size() << endl;

			/* banded HMM align */
			hmm.calcViterbiScores(seqVscore, seqVpath);
			float minCost = hmm.buildViterbiTrace(seqVscore, seqVpath);
			if(minCost == inf) {
				cerr << "Unable to align read " << read.getId() << "No Viterbi path found" << endl;
				continue;
			}
			if(!hmm.isValidAlignPath(seqVpath)) {
				cerr << "Invalid align path: length: " << seqVpath.alnPath.length() << " path: " << seqVpath.alnPath << endl
					 << "alnStart: " << seqVpath.alnStart << " alnEnd: " << seqVpath.alnEnd << " alnLen: "
					 << " alnLen: " << seqVpath.alnEnd - seqVpath.alnStart + 1
					 << " alnFrom: " << seqVpath.alnFrom << " alnTo: " << seqVpath.alnTo << endl;
				abort();
			}
			cerr << "maxScore: " << minCost << endl;
			cerr << "alnStart: " << seqVpath.alnStart << endl;
			cerr << "alnEnd: " << seqVpath.alnEnd << endl;

			/* find seqStart and seqEnd */
			int csStart = hmm.getCSLoc(seqVpath.alnStart) - 1;
			int csEnd = hmm.getCSLoc(seqVpath.alnEnd) - 1;
			cerr << "csStart: " << csStart << " csEnd: " << csEnd << " csLen: " << (csEnd - csStart + 1) << endl;

			PrimarySeq aln = hmm.buildGlobalAlignSeq(seqVscore, seqVpath);
			DigitalSeq seq(aln);
			if(seq.length() != ptu.numAlignSites()) {
				cerr << "aln.length: " << aln.length() << " seq.length: " << seq.length() << " numSite: " << ptu.numAlignSites() << endl;
				continue;
			}

			/* phylogenetic placement */
			double maxLoglik = EGriceLab::infV;
			PTUnrooted::PTUNodePtr bestNode;
			string bestAnnotation;
			/* Use a SLA (leaf-ancestor) search algorithm */
			const vector<PTUnrooted::PTUNodePtr>& leafHits = ptu.getLeafHits(getLeafHitsById(seqIdHits, id2leaf), seq, maxDist, csStart, csEnd);
			infoLog << "Found " << leafHits.size() << " leaf nodes for " << read.getId() << endl;
			if(leafHits.empty())
				continue;

			NodeSet nodeSeen;
			for(vector<PTUnrooted::PTUNodePtr>::const_iterator leaf = leafHits.begin(); leaf != leafHits.end(); ++leaf) {
				PTUnrooted::PTUNodePtr node = *leaf;
				double prevlik = EGriceLab::infV;
				while(!node->isRoot() && nodeSeen.find(node) == nodeSeen.end()) {
					/* place the read here */
					PTUnrooted subtree = ptu.copySubTree(node, node->getParent());
					const PTUnrooted::PTUNodePtr& v = subtree.getNode(0);
					const PTUnrooted::PTUNodePtr& u = subtree.getNode(1);
					subtree.placeSeq(seq, u, v, csStart, csEnd);
					const PTUnrooted::PTUNodePtr& r = subtree.getNode(2);
					const PTUnrooted::PTUNodePtr& n = subtree.getNode(3);
					double wrn = subtree.getBranchLength(r, n);
					double treeLik = subtree.treeLoglik(csStart, csEnd);
					nodeSeen.insert(node);

					if(treeLik > maxLoglik) {
						maxLoglik = treeLik;
						bestNode = node;
						bestAnnotation = n->getAnnotation();
					}
					if(treeLik <= prevlik)
						break;

					node = node->getParent();
					prevlik = treeLik;
				} /* end of this cascade */
			} /* end each leaf */
			out << read.getId() << "\t" << read.getDesc() << "\t"
				<< bestNode->getParent()->getId() << "->" << bestNode->getId() << "\t"
				<< bestAnnotation << "\t"
				<< csStart << "-" << csEnd << "\t" << endl;
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
