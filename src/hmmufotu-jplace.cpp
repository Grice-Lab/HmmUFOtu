/*
 * hmmufotu-jplace.cpp
 *  generate JPlace placement summary from hmmufotu assignment files
 *
 *  Created on: Jan 30, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <cctype>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <limits>
#include <vector>
#include <map>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp> /* for boost string join */
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp> /* basic boost streams */
#include <boost/iostreams/device/file.hpp> /* file sink and source */
#include <boost/iostreams/filter/zlib.hpp> /* for zlib support */
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp> /* for bzip2 support */
#include <json/json.h> /* jsoncpp support */
#include "HmmUFOtu.h"
#include "HmmUFOtu_main.h"

using namespace std;
using namespace EGriceLab;
using namespace EGriceLab::HmmUFOtu;
using namespace Eigen;

/* default values */
static const double DEFAULT_MIN_Q = 0;
static const double DEFAULT_MIN_ALN_IDENTITY = 0;
static const double DEFAULT_MIN_HMM_IDENTITY = 0;
static const int JPLACE_VERSION = 3;
static const char *field_names[] = { "edge_num", "likelihood", "like_weight_ratio", "distal_length", "proximal_length", "pendant_length" };
static const string TREE_NODE_NAME = "tree";
static const string PLACEMENT_LIST_NODE_NAME = "placements";
static const string PLACEMENT_NODE_NAME = "p";
static const string READNAME_NODE_NAME = "n";
static const string VERSION_NODE_NAME = "version";
static const string FIELD_NODE_NAME = "fields";
static const string INVOCATION_NODE_NAME = "invocation";
static const string METADATA_NODE_NAME = "metadata";
static const string SM_NODE_NAME = "substitution_model";
static const string VAR_NODE_NAME = "among_site_rate_variation";
static const string ANNO_NODE_NAME = "node_taxonomy_annotations";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Generate JPlace (JSON phylogenetic-placement) file from HmmUFOtu taxonomy assignment files" << endl
		 <<	"Warning:  this program need go store all records in memory and could take large amount of RAM, proceed with causion" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed file";
	#endif
	cerr << "Usage:    " << progName << "  <HmmUFOtu-DB> <(INFILE [INFILE2 ...]> [options]" << endl
		 << "INFILE          INFILE         : assignment file(s) from hmmufotu" << ZLIB_SUPPORT << endl
		 << "Options:    -o  FILE           : write the jplace output to FILE instead of stdout" << endl
		 << "            -q  DBL            : minimum qPlace score (negative log10 posterior error rate) required [" << DEFAULT_MIN_Q << "]" << endl
		 << "            --aln-iden  DBL    : minimum alignment identity required for assignment result [" << DEFAULT_MIN_ALN_IDENTITY << "]" << endl
		 << "            --hmm-iden  DBL    : minimum profile-HMM identity required for assignment result [" << DEFAULT_MIN_HMM_IDENTITY << "]" << endl
		 << "            -sm  FLAG          : report the DNA Substitution Model name in metadata used for the phylogenetic placement" << endl
		 << "            -V|--var  FLAG        : report the among site rate variation model in metadata used for the phylogenetic placement" << endl
		 << "            -a|--anno  FLAG    : report all node taxonomic annotations in metadata" << endl
		 << "            -v  FLAG           : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version          : show program version and exit" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, hmmFn, ptuFn;
	vector<string> inFiles;
	string outFn;
	ifstream hmmIn, ptuIn;
	ofstream of;

	Json::Value jptree; /* create the root */

	double minQ = DEFAULT_MIN_Q;
	double minAlnIden = DEFAULT_MIN_ALN_IDENTITY;
	double minHmmIden = DEFAULT_MIN_HMM_IDENTITY;
	bool showSm = false;
	bool showVar = false;
	bool showAnno = false;

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

	if(!(cmdOpts.numMainOpts() > 1)) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	dbName = cmdOpts.getMainOpt(0);
	for(int i = 1; i < cmdOpts.numMainOpts(); ++i) {
		string fn = cmdOpts.getMainOpt(i);
		inFiles.push_back(fn);
	}

	if(cmdOpts.hasOpt("-o"))
		outFn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("-q"))
		minQ = ::atof(cmdOpts.getOptStr("-q"));

	if(cmdOpts.hasOpt("--aln-iden"))
		minAlnIden = ::atof(cmdOpts.getOptStr("--aln-iden"));
	if(cmdOpts.hasOpt("--hmm-iden"))
		minHmmIden = ::atof(cmdOpts.getOptStr("--hmm-iden"));

	if(cmdOpts.hasOpt("-sm"))
		showSm = true;
	if(cmdOpts.hasOpt("-V") || cmdOpts.hasOpt("--var"))
		showVar = true;

	if(cmdOpts.hasOpt("-a") || cmdOpts.hasOpt("--anno"))
		showAnno = true;

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* validate options */
	if(!(minQ >= 0)) {
		cerr << "-q must be non-negative" << endl;
		return EXIT_FAILURE;
	}

	/* set filenames */
	hmmFn = dbName + HMM_FILE_SUFFIX;
	ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;

	/* open inputs */
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

	/* load database */
	BandedHMMP7 hmm;
	hmmIn >> hmm;
	if(hmmIn.bad()) {
		cerr << "Unable to read HMM profile '" << hmmFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "HMM profile read" << endl;

	if(loadProgInfo(ptuIn).bad())
		return EXIT_FAILURE;
	PTUnrooted ptu;
	ptu.load(ptuIn); /* only load the tree topology, ignore loglik for space */
	if(ptuIn.bad()) {
		cerr << "Unable to load Phylogenetic tree data '" << ptuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	infoLog << "Phylogenetic tree loaded" << endl;

	/* add tree structure */
	jptree[TREE_NODE_NAME] = ptu.toJPlaceTreeStr(ptu.getRoot()) + ";";
	Json::Value placements_list;

	/* process input files */
	for(vector<string>::const_iterator infn = inFiles.begin(); infn != inFiles.end(); ++infn) {
		infoLog << "Processing " << *infn << " ..." << endl;
		boost::iostreams::filtering_istream in;

#ifdef HAVE_LIBZ
		if(StringUtils::endsWith(*infn, GZIP_FILE_SUFFIX))
			in.push(boost::iostreams::gzip_decompressor());
		else if(StringUtils::endsWith(*infn, BZIP2_FILE_SUFFIX))
			in.push(boost::iostreams::bzip2_decompressor());
		else { }
#endif

		in.push(boost::iostreams::file_source(*infn));
		if(in.bad()) {
			cerr << "Unable to open assignment input file '" << *infn << "' " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}

		/* check program info */
		if(readProgInfo(in).bad())
			return EXIT_FAILURE;

		TSVScanner tsvIn(in, true);
		while(tsvIn.hasNext()) {
			const TSVRecord& record = tsvIn.nextRecord();

			string rid = record.getFieldByName("id");
			int csStart = ::atoi(record.getFieldByName("CS_start").c_str());
			int csEnd = ::atoi(record.getFieldByName("CS_end").c_str());
			const string& aln = record.getFieldByName("alignment");
			const string& branch_id = record.getFieldByName("branch_id");
			double branch_ratio = ::atof(record.getFieldByName("branch_ratio").c_str());
			const long taxon_id = ::atol(record.getFieldByName("taxon_id").c_str());
			double annoDist = ::atof(record.getFieldByName("anno_dist").c_str());
			double loglik = ::atof(record.getFieldByName("loglik").c_str());
			double q = ::atof(record.getFieldByName("Q_placement").c_str());

			if(taxon_id >= 0 && q >= minQ
					&& alignIdentity(AlphabetFactory::nuclAbc, aln, csStart - 1, csEnd -1)
			&& hmmIdentity(hmm, aln, csStart - 1, csEnd - 1)) { /* a valid assignment */
				long pNodeId = 0;
				long cNodeId = 0;
				sscanf(branch_id.c_str(), "%d->%d", &cNodeId, &pNodeId);
				PTUnrooted::PTUNodePtr cNode = ptu.getNode(cNodeId);
				PTUnrooted::PTUNodePtr pNode = ptu.getNode(pNodeId);
				JPlace place(ptu.getEdgeID(cNode, pNode), rid, ptu.getBranchLength(cNode, pNode),
						branch_ratio, loglik, annoDist, q);

				Json::Value place_node;

				/* add a one-row placement matrix */
				Json::Value pmatrix;
				pmatrix[0].append(place.edgeID);
				pmatrix[0].append(place.likelihood);
				pmatrix[0].append(place.like_ratio);
				pmatrix[0].append(place.distal_length);
				pmatrix[0].append(place.proximal_length);
				pmatrix[0].append(place.pendant_length);

				place_node[PLACEMENT_NODE_NAME] = pmatrix;

				/* add a one-element read name list */
				Json::Value read_node;
				read_node.append(place.readName);
				place_node[READNAME_NODE_NAME] = read_node;

				/* add this place_node to place_list */
				placements_list.append(place_node);
			} /* end if */
		} /* end each record */
	} /* end eachi file */
	/* add placements_list */
	jptree[PLACEMENT_LIST_NODE_NAME] = placements_list;

	/* add mendatory metadata */
	jptree[VERSION_NODE_NAME] = JPLACE_VERSION;
	for(const char** name = field_names; name != field_names + sizeof(field_names)/sizeof(*field_names); ++name)
		jptree[FIELD_NODE_NAME].append(*name);

	/* add optional metadata */
	Json::Value metadata;
	metadata[INVOCATION_NODE_NAME] = cmdOpts.getCmdStr();
	if(showSm)
		metadata[SM_NODE_NAME] = ptu.getModel()->modelType();

	if(showSm)
		metadata[VAR_NODE_NAME] = ptu.isVar() ? "Discrete Gamma model" : "none";

	if(showAnno) {
		const vector<PTUnrooted::PTUNodePtr>& allNodes = ptu.getNodes();
		Json::Value anno_list;
		for(vector<PTUnrooted::PTUNodePtr>::const_iterator node = allNodes.begin(); node != allNodes.end(); ++node)
			anno_list[boost::lexical_cast<string>((*node)->getId())] = (*node)->getAnno();
		metadata[ANNO_NODE_NAME] = anno_list;
	}
	jptree[METADATA_NODE_NAME] = metadata;


	/* write jptree */
    out << jptree << endl;
}



