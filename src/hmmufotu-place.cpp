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
static const string TBL_HEADER = "id\tdesc\tbranch_id\tbranch_name\tphylogenetic_annotation";

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Phylogenetic placement based taxonamy assignment for multiple-aligned sequences" << endl
		 << "Usage:    " << progName << "  <HmmUFOtu-DB> <MSA-INFILE> [options]" << endl
		 << "MSA-FILE  FILE                 : multiple-sequence aligned input" << endl
		 << "Options:    -o  FILE           : write output to FILE instead of stdout" << endl
		 << "            -T  FILE           : in addition to the placement output, write the taxa summary table to TAXA-FILE" << endl
		 << "            -d|--pdist  DBL    : maximum p-dist between placed read and observed leaves during the optimal search [" << DEFAULT_MAX_PDIST << "]" << endl
		 << "            -q  DBL            : minimum Q-value (negative log10 of the error) required for a read placement [" << DEFAULT_MIN_Q << "]" << endl
		 << "            -v  FLAG           : enable verbose information" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string infn, seqfn, outfn, taxafn;
	string fmt;
	ifstream ptuIn;
	ofstream of, taxaOut;

	PTUnrooted ptu;
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

	infn = cmdOpts.getMainOpt(0);
	if(!StringUtils::endsWith(infn, MSA_FILE_SUFFIX))
		infn += MSA_FILE_SUFFIX;

	seqfn = cmdOpts.getMainOpt(1);
	if(StringUtils::endsWith(seqfn, ".fasta") || StringUtils::endsWith(seqfn, ".fas")
		|| StringUtils::endsWith(seqfn, ".fa") || StringUtils::endsWith(seqfn, ".fna"))
		fmt = "fasta";
	else if(StringUtils::endsWith(seqfn, ".fastq") || StringUtils::endsWith(seqfn, ".fq"))
		fmt = "fastq";
	else {
		cerr << "Unrecognized format of MSA file '" << seqfn << "'" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-o"))
		outfn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("-d"))
		maxDist = ::atof(cmdOpts.getOptStr("-d"));
	if(cmdOpts.hasOpt("--pdist"))
		maxDist = ::atof(cmdOpts.getOptStr("--pdist"));

	if(cmdOpts.hasOpt("-q"))
		minQ = ::atoi(cmdOpts.getOptStr("-q"));

	if(cmdOpts.hasOpt("-v"))
		ENABLE_INFO();

	/* open files */
	ptuIn.open(infn.c_str(), ios_base::in | ios_base::binary);
	if(!ptuIn.is_open()) {
		cerr << "Unable to open '" << infn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	SeqIO seqIn(seqfn, ALPHABET, fmt);

	if(!outfn.empty()) {
		of.open(outfn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to '" << outfn << "'" << endl;
			return EXIT_FAILURE;
		}
	}
	ostream& out = of.is_open() ? of : cout;

	infoLog << "Loading Phylogenetic tree" << endl;
	ptu.load(ptuIn);
	if(ptuIn.bad()) {
		cerr << "Failed to load PTU data from " << infn << endl;
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
		/* Use a LA (leaf-ancestor) search algorithm */
		vector<PTUnrooted::PTUNodePtr> leafHits;
		for(vector<PTUnrooted::PTUNodePtr>::size_type i = 0; i < ptu.numNodes(); ++i) {
			const PTUnrooted::PTUNodePtr& node = ptu.getNode(i);
			if(node->isLeaf() && DNASubModel::pDist(node->getSeq(), seq, start, end + 1) <= maxDist)
				leafHits.push_back(node);
		}
		infoLog << "Found " << leafHits.size() << " leaf nodes for " << read.getId() << endl;
		if(leafHits.empty())
			continue;

		NodeSet nodeSeen;
		for(vector<PTUnrooted::PTUNodePtr>::const_iterator it = leafHits.begin(); it != leafHits.end(); ++it) {
			PTUnrooted::PTUNodePtr node = *it; /* make a copy */
			while(!node->isRoot() && nodeSeen.find(node) == nodeSeen.end()) {
				/* place the read here */
				PTUnrooted subtree = ptu.copySubTree(node, node->getParent());
				subtree.placeSeq(seq, subtree.getNode(1), subtree.getNode(0), start, end);
				double vn = subtree.getBranchLength(subtree.getNode(subtree.numNodes() - 2), subtree.getNode(subtree.numNodes() - 1));
				double tc = subtree.treeLoglik(start, end);

//				infoLog << "Read " << read.getId() << " placed at node " << node->getId() <<
//						" new branch length: " << vn
//						<< " cost: " << tc << endl;
				if(tc > maxLoglik) {
					maxLoglik = tc;
					bestNode = node;
				}
				node = node->getParent();
			} /* end of this cascade */
		} /* end each leaf */
		out << read.getId() << "\t" << read.getDesc() << "\t"
			<< bestNode->getParent()->getId() << "->" << bestNode->getId() << "\t"
			<< bestNode->getParent()->getName() << "->" << bestNode->getName() << "\t"
			<< "" << endl;
	}
}
