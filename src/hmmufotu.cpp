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
#include <boost/iostreams/filtering_stream.hpp> /* basic boost streams */
#include <boost/iostreams/device/file.hpp> /* file sink and source */
#include <boost/iostreams/filter/zlib.hpp> /* for zlib support */
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp> /* for bzip2 support */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "HmmUFOtu.h"
#include "HmmUFOtu_main.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::HmmUFOtu;
using namespace Eigen;

/* default values */
static const double DEFAULT_MAX_DIFF = inf;
static const size_t DEFAULT_MAX_LOCATION = 50;
static const int DEFAULT_SEED_LEN = 20;
static const int MAX_SEED_LEN = 25;
static const int MIN_SEED_LEN = 15;
static const int DEFAULT_SEED_REGION = 50;
static const double DEFAULT_MAX_PLACE_ERROR = 20;
static const int DEFAULT_NUM_SEGMENT = 2;
static const int MIN_NUM_SEGMENT = 2;
static const int MAX_NUM_SEGMENT = 4;
static const int DEFAULT_NUM_THREADS = 1;
static const string ALIGN_OUT_FMT = "fasta";
static const string DEFAULT_BRANCH_EST_METHOD = "unweighted";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Ultra-fast microbiome amplicon sequencing read taxonomy assignment and OTU picking tool,"
		 << "based on Consensus-Sequence-FM-index (CSFM-index) powered HMM alignment"
		 << " and Seed-Estimate-Place (SEP) local phylogenetic placement" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed file";
	#endif

	cerr << "Usage:    " << progName << "  <HmmUFOtu-DB> <READ-FILE1> [READ-FILE2] [options]" << endl
		 << "READ-FILE1  FILE                 : sequence read file for the assembled/forward read" << ZLIB_SUPPORT << endl
		 << "READ-FILE2  FILE                 : sequence read file for the reverse read" << ZLIB_SUPPORT << endl
		 << "Options:    -o  FILE             : write the assignment output to FILE instead of stdout" << ZLIB_SUPPORT << endl
		 << "            -a  FILE             : in addition to the assignment output, write the read alignment in " << ALIGN_OUT_FMT << " format" << ZLIB_SUPPORT << endl
		 << "            --fmt  STR           : read file format (applied to all read files), supported format: 'fasta', 'fastq'" << endl
		 << "            -L|--seed-len  INT   : seed length used for banded-Hmm search [" << DEFAULT_SEED_LEN << "]" << endl
		 << "            -R  INT              : size of 5'/3' seed region for finding seed matches for CSFM-index [" << DEFAULT_SEED_REGION << "]" << endl
		 << "            -s  FLAG             : assume READ-FILE1 is single-end read instead of assembled read, if no READ-FILE2 provided" << endl
		 << "            -N  INT              : max # of seed nodes used in the 'Seed' stage of SEP algorithm [" << DEFAULT_MAX_LOCATION << "]" << endl
		 << "            -d  DBL              : max p-dist difference allowed for sub-optimal seeds compared the best seed in the 'Estimate' stage of SEP algorithm [" << DEFAULT_MAX_DIFF << "]" << endl
		 << "            -e  DBL              : max placement error for filtering placements between the 'Estimate' and 'Place' stages of the SEP algorithm [" << DEFAULT_MAX_PLACE_ERROR << "]" << endl
		 << "            -m|--method  STR     : branch length estimating method during the estimated-placement stage, must be one of 'unweighted' or 'weighted' [" << DEFAULT_BRANCH_EST_METHOD << "]" << endl
		 << "            --ML  FLAG           : use maximum likelihood in phylogenetic placement, do not calculate posterior p-values, this will ignore -q and --prior options" << endl
		 << "            --prior  STR         : method for calculating prior probability of a placement, either 'uniform' (uniform prior) or 'height' (rooted distance to leaves)" << endl
		 << "            -C|--chimera  FLAG   : enable a chimera sequence checking procedure before the final 'Place' stage in the SEP algorithm using a segment re-estimation method" << endl
		 << "            --num-seg  INT       : number of segments used in chimera checking procedure [" << DEFAULT_NUM_SEGMENT << "]" << endl
		 << "            --chimera-out  FILE  : keep read assignment results of chimera reads in FILE" << ZLIB_SUPPORT << endl
		 << "            -S|--seed  INT       : random seed used for CSFM-index seed searches, for debug only" << endl
#ifdef _OPENMP
		 << "            -p|--process INT     : number of threads/cpus used for parallel processing" << endl
#endif
		 << "            --align-only  FLAG   : only align the read but not try to place it into the tree, this will make " + progName + " behaviors like an HMM aligner" << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	/* filenames */
	string dbName, fwdFn, revFn, msaFn, csfmFn, hmmFn, ptuFn;
	string outFn, alnFn;
	string chiOutFn;
	/* input */
	ifstream msaIn, csfmIn, hmmIn, ptuIn;
	boost::iostreams::filtering_istream fwdIn, revIn;
	/* output */
	boost::iostreams::filtering_ostream out, alnOut;
	boost::iostreams::filtering_ostream chiOut;
	/* other */
	string seqFmt; /* seq file format */
	string estMethod = DEFAULT_BRANCH_EST_METHOD;
	SeqIO fwdSeqI, revSeqI, alnSeqO;

	bool isAssembled = true; /* assume assembled seq if not paired-end */
	bool alignOnly = false;
	BandedHMMP7::align_mode mode;

	int seedLen = DEFAULT_SEED_LEN;
	int seedRegion = DEFAULT_SEED_REGION;
	double maxDiff = DEFAULT_MAX_DIFF;
	int maxLocs = DEFAULT_MAX_LOCATION;
	double maxError = DEFAULT_MAX_PLACE_ERROR;
	bool onlyML = false;
	PTUnrooted::PRIOR_TYPE myPrior = PTUnrooted::UNIFORM;
	bool checkChimera = false;
	int numSeg = DEFAULT_NUM_SEGMENT;

	int nThreads = DEFAULT_NUM_THREADS;

	unsigned seed = time(NULL); // using time as default seed

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

	if(cmdOpts.hasOpt("--fmt"))
		seqFmt = cmdOpts.getOpt("--fmt");

	if(cmdOpts.hasOpt("-L"))
		seedLen = ::atoi(cmdOpts.getOptStr("-L"));
	if(cmdOpts.hasOpt("--seed-len"))
		seedLen = ::atoi(cmdOpts.getOptStr("--seed-len"));

	if(cmdOpts.hasOpt("-R"))
		seedRegion = ::atoi(cmdOpts.getOptStr("-R"));

	if(cmdOpts.hasOpt("-s"))
		isAssembled = false;

	if(cmdOpts.hasOpt("-d"))
		maxDiff = ::atof(cmdOpts.getOptStr("-d"));

	if(cmdOpts.hasOpt("-N"))
		maxLocs = ::atoi(cmdOpts.getOptStr("-N"));

	if(cmdOpts.hasOpt("-e"))
		maxError = ::atof(cmdOpts.getOptStr("-e"));

	if(cmdOpts.hasOpt("-m"))
		estMethod = cmdOpts.getOpt("-m");
	if(cmdOpts.hasOpt("--method"))
		estMethod = cmdOpts.getOpt("--method");

	if(cmdOpts.hasOpt("--ML"))
		onlyML = true;

	if(cmdOpts.hasOpt("--prior")) {
		if(cmdOpts.getOpt("--prior") == "uniform")
			myPrior = PTUnrooted::UNIFORM;
		else if(cmdOpts.getOpt("--prior") == "height")
			myPrior = PTUnrooted::HEIGHT;
		else {
			cerr << "Unsupported prior specified, check the --prior option" << endl;
			return EXIT_FAILURE;
		}
	}

	if(cmdOpts.hasOpt("-C") || cmdOpts.hasOpt("--chimera")) {
		checkChimera = true;
		if(cmdOpts.hasOpt("--num-seg"))
			numSeg = ::atof(cmdOpts.getOptStr("--num-seg"));
		if(cmdOpts.hasOpt("--chimera-out"))
			chiOutFn = cmdOpts.getOpt("--chimera-out");
	}

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

	/* guess fwdSeq format */
	if(seqFmt.empty()) {
		string seqPre = fwdFn;
		StringUtils::removeEnd(seqPre, GZIP_FILE_SUFFIX);
		StringUtils::removeEnd(seqPre, BZIP2_FILE_SUFFIX);
		seqFmt = SeqUtils::guessSeqFileFormat(seqPre);
	}
	if(!(seqFmt == "fasta" || seqFmt == "fastq")) {
		cerr << "Unsupported sequence format '" << seqFmt << "'" << endl;
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
	if(!(maxDiff >= 0)) {
		cerr << "-d must be non-negative" << endl;
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
	if(!(MIN_NUM_SEGMENT <= numSeg && numSeg <= MAX_NUM_SEGMENT)) {
		cerr << "--num-seg must be in [" << MIN_NUM_SEGMENT << ", " << MAX_NUM_SEGMENT << "]" << endl;
		return EXIT_FAILURE;
	}

#ifdef _OPENMP
	if(!(nThreads > 0)) {
		cerr << "-p|--process must be positive" << endl;
		return EXIT_FAILURE;
	}
	omp_set_num_threads(nThreads);
#endif

	bool isSingle = revFn.empty();
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

#ifdef HAVE_LIBZ
	if(StringUtils::endsWith(fwdFn, GZIP_FILE_SUFFIX))
		fwdIn.push(boost::iostreams::gzip_decompressor());
	else if(StringUtils::endsWith(fwdFn, BZIP2_FILE_SUFFIX))
		fwdIn.push(boost::iostreams::bzip2_decompressor());
	else { }
#endif
	if(fwdIn.empty()) /* not zipped */
		fwdIn.push(boost::iostreams::file_source(fwdFn));
	else
		fwdIn.push(boost::iostreams::file_source(fwdFn, std::ios_base::in | std::ios_base::binary));
	if(fwdIn.bad()) {
		cerr << "Unable to open forward seq file '" << fwdFn << "' " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(!revFn.empty()) {
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(revFn, GZIP_FILE_SUFFIX))
			revIn.push(boost::iostreams::gzip_decompressor());
		else if(StringUtils::endsWith(revFn, BZIP2_FILE_SUFFIX))
			revIn.push(boost::iostreams::bzip2_decompressor());
		else { }
#endif
		revIn.push(boost::iostreams::file_source(revFn));
		if(revIn.bad()) {
			cerr << "Unable to open reverse seq file '" << revFn << "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* open outputs */
#ifdef HAVE_LIBZ
	if(StringUtils::endsWith(outFn, GZIP_FILE_SUFFIX)) /* empty outFn won't match */
		out.push(boost::iostreams::gzip_compressor());
	else if(StringUtils::endsWith(outFn, BZIP2_FILE_SUFFIX)) /* empty outFn won't match */
		out.push(boost::iostreams::bzip2_compressor());
	else { }
#endif
	if(!outFn.empty())
		out.push(boost::iostreams::file_sink(outFn));
	else
		out.push(std::cout);
	if(out.bad()) {
		cerr << "Unable to write to "
				<< (!outFn.empty() ? " out file '" + outFn + "' " : "stdout ")
				<< ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(!alnFn.empty()) {
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(alnFn, GZIP_FILE_SUFFIX))
			alnOut.push(boost::iostreams::gzip_compressor());
		else if(StringUtils::endsWith(alnFn, BZIP2_FILE_SUFFIX))
			alnOut.push(boost::iostreams::bzip2_compressor());
		else { }
#endif
		alnOut.push(boost::iostreams::file_sink(alnFn));
		if(alnOut.bad()) {
			cerr << "Unable to write to align file '" << alnFn << "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	if(!chiOutFn.empty()) {
#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(chiOutFn, GZIP_FILE_SUFFIX))
			chiOut.push(boost::iostreams::gzip_compressor());
		else if(StringUtils::endsWith(outFn, BZIP2_FILE_SUFFIX))
			chiOut.push(boost::iostreams::bzip2_compressor());
		else { }
#endif
		chiOut.push(boost::iostreams::file_sink(chiOutFn));
		if(chiOut.bad()) {
			cerr << "Unable to write to '" + chiOutFn + "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	/* loading database files */
	if(loadProgInfo(msaIn).bad())
		return EXIT_FAILURE;
	MSA msa;
	msa.load(msaIn);
	if(msaIn.bad()) {
		cerr << "Failed to load MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	int csLen = msa.getCSLen();
	infoLog << "MSA loaded" << endl;

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
	const DegenAlphabet* abc = hmm.getNuclAbc();

	/* prepare SeqIO */
	fwdSeqI.reset(dynamic_cast<istream*> (&fwdIn), abc, seqFmt);
	if(!revFn.empty())
		revSeqI.reset(dynamic_cast<istream*> (&revIn), abc, seqFmt);

	if(!alnFn.empty())
		alnSeqO.reset(dynamic_cast<ostream*> (&alnOut), abc, ALIGN_OUT_FMT);

	debugLog << "Sequence input and output prepared" << endl;

	if(loadProgInfo(csfmIn).bad())
		return EXIT_FAILURE;
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

	if(loadProgInfo(ptuIn).bad())
		return EXIT_FAILURE;
	PTUnrooted ptu;
	if(!alignOnly) {
		ptu.load(ptuIn);
		if(ptuIn.bad()) {
			cerr << "Unable to load Phylogenetic tree data '" << ptuFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		infoLog << "Phylogenetic tree loaded" << endl;
	}

	/* configure HMM mode */
	hmm.setSequenceMode(mode);
	hmm.wingRetract();

	infoLog << "Processing read ..." << endl;
	/* process reads and output */
	writeProgInfo(out, string(" taxonomy assignment generated by ") + argv[0]);
	out << "# command: "<< cmdOpts.getCmdStr() << endl;
	out << "id\tdescription\t" << BandedHMMP7::HmmAlignment::TSV_HEADER << "\t" + PTUnrooted::PTPlacement::TSV_HEADER << endl;
	if(chiOut.is_complete())
		chiOut << "id\tdescription\t" << BandedHMMP7::HmmAlignment::TSV_HEADER << "\t" + PTUnrooted::PTPlacement::TSV_HEADER << endl;

#pragma omp parallel
	{
#pragma omp single
		{
			while(fwdSeqI.hasNext() && (revFn.empty() || revSeqI.hasNext())) {
				string id;
				string desc;
				PrimarySeq fwdRead, revRead;
				bool isPaired = true;
				bool isOriented = true;

				fwdRead = fwdSeqI.nextSeq();
				id = fwdRead.getId();
				desc = fwdRead.getDesc();
				if(!revFn.empty()) { /* paired-ended */
					revRead = revSeqI.nextSeq().revcom();
					assert(fwdRead.getId() == revRead.getId());
				}
#pragma omp task
				{
					BandedHMMP7::HmmAlignment aln;
					/* align fwdRead */
					aln = alignSeq(hmm, csfm, fwdRead, seedLen, seedRegion, mode);
					assert(aln.isValid());
					//						infoLog << "fwd seq aligned: csStart: " << csStart << " csEnd: " << csEnd << " aln: " << aln << endl;
					if(!revFn.empty()) { /* align revRead */
						//							cerr << "Aligning mate: " << revRead.getId() << endl;
						BandedHMMP7::HmmAlignment revAln = alignSeq(hmm, csfm, revRead, seedLen, seedRegion, mode);
						assert(revAln.isValid());
						//							infoLog << "rev seq aligned: revStart: " << revStart << " revEnd: " << revEnd << " aln: " << revAln << endl;
						if(!(aln.csStart <= revAln.csStart && aln.csEnd <= revAln.csEnd)) /* bad paired-end alignment */
							isOriented = false;
						else
							aln.merge(revAln); /* merge alignment */
					}

					if(isOriented) { /* a valid alignment */
						/* write the alignment seq to output */
						if(!alnFn.empty()) {
							string desc = fwdRead.getDesc();
							desc += ";csStart=" + boost::lexical_cast<string>(aln.csStart) +
									";csEnd=" + boost::lexical_cast<string>(aln.csEnd) + ";";
#pragma omp critical(writeAln)
							alnSeqO.writeSeq(PrimarySeq(abc, id, aln.align, desc));
						}

						PTUnrooted::PTPlacement bestPlace;

						if(!alignOnly) {
							DigitalSeq seq(abc, id, aln.align);
							/* place seq with seed-estimate-place (SEP) algorithm */
							/* get potential locs */
							vector<PTUnrooted::PTLoc> locs = getSeed(ptu, seq, aln.csStart - 1, aln.csEnd - 1, maxDiff);
							if(locs.size() > maxLocs)
								locs.erase(locs.end() - (locs.size() - maxLocs), locs.end()); /* remove last maxLocs elements */
							//	cerr << "Found " << locs.size() << " potential placement locations" << endl;
							/* estimate placements */
							vector<PTUnrooted::PTPlacement> places = estimateSeq(ptu, seq, locs, estMethod);
							/* filter placements */
							filterPlacements(places, maxError);
							/* accurate placements */
							placeSeq(ptu, seq, places);
							if(onlyML) { /* don't calculate q-values */
								std::sort(places.rbegin(), places.rend(), compareByLoglik); /* sort places decently by real loglik */
							}
							else { /* calculate q-values */
								calcQValues(places, myPrior);
								std::sort(places.rbegin(), places.rend(), compareByQPlace); /* sort places decently by posterior placement probability */
							}

							bestPlace = places[0];

						} /* end if alignOnly */
						/* write main output */
#pragma omp critical(writeAssign)
						out << id << "\t" << desc << "\t"
								<< aln << "\t"
								<< bestPlace << endl;
						if(chiOut.is_complete())
#pragma omp critical(writeChiAssign)
							chiOut << id << "\t" << desc << "\t"
							<< aln << "\t"
							<< bestPlace << endl;
					} /* end valid alignment */
				} /* end task */
			} /* end each read/pair */
		} /* end single */
#pragma omp taskwait
	} /* end parallel */
	/* release resources */
}

