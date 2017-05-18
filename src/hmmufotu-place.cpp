/*
 * hmmufotu-place.cpp
 *  Place aligned 16S RNA reads into a fixed phylogenetic tree pre-built by hmmufotu-build
 *
 *  Created on: Jan 12, 2017
 *      Author: zhengqi
 */

#include <iostream>
#include <fstream>
#include <cfloat>
#include <boost/unordered_set.hpp>
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_phylo.h"

using namespace std;
using namespace EGriceLab;

/** default values */
static const string ALPHABET = "dna";

static const double DEFAULT_MAX_PDIST = 0.03;
static const int DEFAULT_MIN_Q = 0;
static const string TBL_HEADER = "id\tdesc\tbranch_id\tbranch_name\tplace_region\tphylogenetic_annotation";

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Phylogenetic placement based taxonamy assignment for multiple-aligned sequences" << endl
		 << "Usage:    " << progName << "  <HmmUFOtu-DB> <MSA-INFILE> [options]" << endl
		 << "MSA-FILE  FILE                 : multiple-sequence aligned input" << endl
		 << "Options:    -o  FILE           : write output to FILE instead of stdout" << endl
		 << "            -T  FILE           : in addition to the placement output, write the taxon summary table to TAXON-FILE" << endl
		 << "            -d|--pdist  DBL    : maximum p-dist between placed read and observed leaves during the optimal search [" << DEFAULT_MAX_PDIST << "]" << endl
		 << "            -q  DBL            : minimum Q-value (negative log10 of the error) required for a read placement [" << DEFAULT_MIN_Q << "]" << endl
		 << "            -v  FLAG           : enable verbose information" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string inFn, seqFn, outFn, taxonFn;
	string fmt;
	ifstream ptuIn;
	ofstream of, taxonOut;

	double maxDist = DEFAULT_MAX_PDIST;
	double minQ = DEFAULT_MIN_Q;

	typedef boost::unordered_set<PTUnrooted::PTUNodePtr> NodeSet;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(cmdOpts.numMainOpts() != 2) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	inFn = cmdOpts.getMainOpt(0);
	if(!StringUtils::endsWith(inFn, PHYLOTREE_FILE_SUFFIX))
		inFn += PHYLOTREE_FILE_SUFFIX;

	seqFn = cmdOpts.getMainOpt(1);
	if(StringUtils::endsWith(seqFn, ".fasta") || StringUtils::endsWith(seqFn, ".fas")
		|| StringUtils::endsWith(seqFn, ".fa") || StringUtils::endsWith(seqFn, ".fna"))
		fmt = "fasta";
	else if(StringUtils::endsWith(seqFn, ".fastq") || StringUtils::endsWith(seqFn, ".fq"))
		fmt = "fastq";
	else {
		cerr << "Unrecognized format of MSA file '" << seqFn << "'" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-o"))
		outFn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("-d"))
		maxDist = ::atof(cmdOpts.getOptStr("-d"));
	if(cmdOpts.hasOpt("--pdist"))
		maxDist = ::atof(cmdOpts.getOptStr("--pdist"));

	if(cmdOpts.hasOpt("-q"))
		minQ = ::atoi(cmdOpts.getOptStr("-q"));
	if(minQ < 0) {
		cerr << "-q must be non-negative" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* open files */
	ptuIn.open(inFn.c_str(), ios_base::in | ios_base::binary);
	if(!ptuIn.is_open()) {
		cerr << "Unable to open '" << inFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	SeqIO seqIn(seqFn, ALPHABET, fmt);

	if(!outFn.empty()) {
		of.open(outFn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to '" << outFn << "'" << endl;
			return EXIT_FAILURE;
		}
	}
	ostream& out = of.is_open() ? of : cout;

	infoLog << "Loading Phylogenetic tree" << endl;
	PTUnrooted ptu;
	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Failed to load PTU data from " << inFn << endl;
		return EXIT_FAILURE;
	}
	else
		infoLog << "Phylogenetic tree data loaded, numNode: " << ptu.numNodes() << " numSites:" << ptu.numAlignSites() << endl;

	const int csLen = ptu.numAlignSites();

	out << TBL_HEADER << endl;

	/* place each seq */
	while(seqIn.hasNext()) {
		const PrimarySeq& read = seqIn.nextSeq();
		if(read.length() != csLen) {
			warningLog << "Warning: aligned read " << read.getId() << " with a length " << read.length() << " not matching the Phylogenetic tree" << endl;
			continue;
		}
		DigitalSeq seq(read);

		/* find the evaluation region */
		int start, end;
		for(start = 0; start < seq.length(); ++start)
			if(seq[start] >= 0) /* valid symbol */
				break;

		for(end = seq.length() - 1; end >= 0; --end)
			if(seq[end] >= 0) /* valid symbol */
				break;

		if(!(start <= end)) {
			warningLog << "Warning: empty read " << read.getId() << " with all gaps found, ignore" << endl;
			continue;
		}

		double maxLoglik = EGriceLab::infV;
		PTUnrooted::PTUNodePtr bestNode;
		string bestAnnotation;
		/* Use a LA (leaf-ancestor) search algorithm */
		vector<PTUnrooted::PTUNodePtr> leafHits;
		for(vector<PTUnrooted::PTUNodePtr>::size_type i = 0; i < ptu.numNodes(); ++i) {
			const PTUnrooted::PTUNodePtr& node = ptu.getNode(i);
			if(node->isLeaf() && SeqUtils::pDist(node->getSeq(), seq, start, end + 1) <= maxDist)
				leafHits.push_back(node);
		}
		infoLog << "Found " << leafHits.size() << " leaf nodes for " << read.getId() << endl;
		if(leafHits.empty())
			continue;

		NodeSet nodeSeen;
		for(vector<PTUnrooted::PTUNodePtr>::const_iterator it = leafHits.begin(); it != leafHits.end(); ++it) {
			PTUnrooted::PTUNodePtr node = *it;
			double prevlik = EGriceLab::infV;
			while(!node->isRoot() && nodeSeen.find(node) == nodeSeen.end()) {
				/* place the read here */
				PTUnrooted subtree = ptu.copySubTree(node, node->getParent());
				const PTUnrooted::PTUNodePtr& v = subtree.getNode(0);
				const PTUnrooted::PTUNodePtr& u = subtree.getNode(1);
				subtree.placeSeq(seq, u, v, start, end);
				const PTUnrooted::PTUNodePtr& r = subtree.getNode(2);
				const PTUnrooted::PTUNodePtr& n = subtree.getNode(3);
				double wrn = subtree.getBranchLength(r, n);
				double treeLik = subtree.treeLoglik(start, end);

				nodeSeen.insert(node);
				double dbest = treeLik - maxLoglik;
				if(treeLik > maxLoglik) {
					maxLoglik = treeLik;
					bestNode = node;
					bestAnnotation = n->getTaxon();
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
			<< start << "-" << end << "\t" << endl;
	}
}
