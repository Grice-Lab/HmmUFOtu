/*
 ============================================================================
 Name        : hmmufotu
 Author      : Qi Zheng
 Version     : v1.1
 Description : Main program of the HmmUFOtu project
 ============================================================================
 */

#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <boost/algorithm/string.hpp> /* for boost string split and join */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "HmmUFOtu.h"
#include "HmmUFOtu_main.h"

using namespace std;
using namespace EGriceLab;
using namespace Eigen;

/* default values */
static const double DEFAULT_MAX_PDIST = inf;
static const size_t DEFAULT_MAX_LOCATION = 50;
static const int DEFAULT_MIN_Q = 0;
static const int DEFAULT_SEED_LEN = 20;
static const int MAX_SEED_LEN = 25;
static const int MIN_SEED_LEN = 15;
static const int DEFAULT_SEED_REGION = 50;
static const double DEFAULT_MAX_PLACE_ERROR = 20;
static const int DEFAULT_NUM_THREADS = 1;

static const string PLACE_HEADER = "id\tdescription\tCS_start\tCS_end\talignment\talign_identity\thmm_identity\tbranch_id\tbranch_ratio\ttaxon_id\ttaxon_anno\tanno_dist\tloglik\tQ_placement\tQ_taxon";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Ultra-fast 16S read OTU assignment using profile-HMM and phylogenetic placement" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <HmmUFOtu-DB> <READ-FILE1> [READ-FILE2] [options]" << endl
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
		 << "            --align-only  FLAG : only align the read but not try to place it into the tree, this will make " + progName + " behaviors like a HMM aligner" << endl
		 << "            -v  FLAG           : enable verbose information, you may set multiple -v for more details" << endl
		 << "            -h|--help          : print this message and exit" << endl;
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
	bool alignOnly = false;
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
		printIntro();
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

	if(cmdOpts.hasOpt("--align-only"))
		alignOnly = true;
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
	if(!alignOnly) {
		ptu.load(ptuIn);
		if(ptuIn.bad()) {
			cerr << "Unable to load Phylogenetic tree data '" << ptuFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		infoLog << "Phylogenetic tree loaded" << endl;
	}

	const DegenAlphabet* abc = hmm.getNuclAbc();

	/* configure HMM mode */
	hmm.setSequenceMode(mode);
	hmm.wingRetract();

	/* process reads and output */
	out << "# Taxonomy assignment generated by: " << progName << " " << progVersion << endl;
	out << "# command: "<< cmdOpts.getCmdStr() << endl;
	out << PLACE_HEADER << endl;
#pragma omp parallel shared(csfm, hmm, ptu, fwdIn, revIn)
#pragma omp single
	{
		while(fwdIn.hasNext() && (!revIn.is_open() || revIn.hasNext())) {
			string id;
			string desc;
			PrimarySeq fwdRead, revRead;
			bool isPaired = true;
			bool isOriented = true;

			fwdRead = fwdIn.nextSeq();
			id = fwdRead.getId();
			desc = fwdRead.getDesc();
			if(revIn.is_open()) { /* paired-ended */
				revRead = revIn.nextSeq().revcom();
				if(fwdRead.getId() != revRead.getId())
					isPaired = false;
			}
			if(!isPaired)
				warningLog << "Warning: Paired-end read doesn't match at " << fwdRead.getId() << " vs. " << revRead.getId() << " , ignored" << endl;
			else
			{
#pragma omp task
				{
					int csStart = 0; /* 1-based consensus start */
					int csEnd = 0;   /* 1-based consensus end */
					string aln;

					/* align fwdRead */
					aln = alignSeq(hmm, csfm, fwdRead, seedLen, seedRegion, mode, csStart, csEnd);
					assert(aln.length() == csLen);
					if(revIn.is_open()) { /* align revRead */
						// cerr << "Aligning mate: " << fwdRead.getId() << endl;
						int revStart = 0;
						int revEnd = 0;
						string revAln = alignSeq(hmm, csfm, revRead, seedLen, seedRegion, mode, revStart, revEnd);
						assert(revAln.length() == csLen);
						if(!(csStart <= revStart && csEnd <= revEnd)) {
							isOriented = false;
#pragma omp critical(warning)
							warningLog << "Warning: Potential incorrectly oriented read found at " << id << " , ignored" << endl;
						}

						/* update alignment */
						csStart = ::min(csStart, revStart);
						csEnd = ::max(csEnd, revEnd);
						BandedHMMP7::mergeWith(aln, revAln);
					}
					assert(1 <= csStart && csStart <= csEnd && csEnd <= csLen);

					if(alnOut.is_open() && isOriented) { /* write the alignment seq to output */
						string desc = fwdRead.getDesc();
						desc += ";csStart=" + boost::lexical_cast<string>(csStart) + ";csEnd=" + boost::lexical_cast<string>(csEnd) + ";";
#pragma omp critical(writeAln)
						alnOut.writeSeq(PrimarySeq(abc, id, aln, fwdRead.getDesc()));
					}
					/* get align and hmm identities */
					double alignIden = alignIdentity(abc, aln, csStart - 1, csEnd - 1);
					double hmmIden = hmmIdentity(hmm, aln, csStart - 1, csEnd - 1);

					PTPlacement bestPlace;
					if(!alignOnly && isOriented) {
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
							std::sort(places.rbegin(), places.rend(), EGriceLab::compareByLoglik); /* sort places decently by estimated loglik */
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
								std::sort(places.rbegin(), places.rend(), EGriceLab::compareByLoglik); /* sort places decently by real loglik */
							}
							else { /* calculate q-values */
								PTPlacement::calcQValues(places, myPrior);
								std::sort(places.rbegin(), places.rend(), EGriceLab::compareByQTaxon); /* sort places decently by QTaxon */
							}

							bestPlace = places.front();
						}
					} /* end if bestOnly */
					/* write main output */
					if(isOriented)
#pragma omp critical(writePlace)
						out << id << "\t" << desc << "\t" << csStart << "\t" << csEnd << "\t"
						<< aln << "\t" << alignIden << "\t" << hmmIden << "\t"
						<< bestPlace << endl;
				} /* end task */
			} /* end else */
		} /* end each read/pair */
#pragma omp taskwait
	} /* end single */
}


