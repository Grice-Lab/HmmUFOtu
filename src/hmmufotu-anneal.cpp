/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/*
 * hmmufotu-anneal.cpp
 *
 *  Created on: May 10, 2017
 *      Author: zhengqi
 *      Version: v1.1
 *      Description : Anneal primer sequences to HmmUFOtu database
 */

#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#include "HmmUFOtu.h"

using namespace std;
using namespace EGriceLab;
using namespace Eigen;

/* default values */
static const double DEFAULT_MIN_IDENTITY = 0.9;
static const int DEFAULT_STRAND = 3;

static const char* TAXON_NAMES[7] = { "KINDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES" };

enum TaxonLevel { KINDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES };

static const string ANNEAL_HEADER = "id\tdescription\tsequence\tstrand\tCS_start\tCS_end\talignment\ttotal_nodes\ttotal_leaves\thit_nodes\thit_leaves\tefficiency_nodes\tefficiency_leaves";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Anneal primer sequences to an HmmUFOtu database and evaluate primer efficiency" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <HmmUFOtu-DB> <SEQ-FILE> [options]" << endl
		 << "SEQ-FILE  FILE                 : primer sequence read file in fasta format, degenerated bases are allowed" << endl
		 << "Options:    -o  FILE           : write the PLACEMENT output to FILE instead of stdout" << endl
		 << "            -i|--identity  DBL : minimum identity between aligned primer sequence and an OTU sequence considered as a good hit [" << DEFAULT_MIN_IDENTITY << "]" << endl
		 << "            -s|--strand  INT   : strand orientation for primers, 1 for forward, 2 for reverse, 3 for auto-detect by best alignment [" << DEFAULT_STRAND << "]" << endl
		 << "            -v  FLAG           : enable verbose information, you may set multiple -v for more details" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, seqFn, msaFn, csfmFn, hmmFn, ptuFn;
	string outFn;
	ifstream msaIn, csfmIn, hmmIn, ptuIn;
	const string seqFmt = "fasta";
	ofstream of;
	ifstream seqIn;
	SeqIO seqI;
	const BandedHMMP7::align_mode mode = BandedHMMP7::GLOBAL;

	double maxDist = 1 - DEFAULT_MIN_IDENTITY;
	int searchStrand = DEFAULT_STRAND;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printIntro();
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}

	if(cmdOpts.numMainOpts() != 2) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}
	dbName = cmdOpts.getMainOpt(0);
	seqFn = cmdOpts.getMainOpt(1);

	if(cmdOpts.hasOpt("-o"))
		outFn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("-i"))
		maxDist = 1 - ::atof(cmdOpts.getOptStr("-i"));
	if(cmdOpts.hasOpt("--identity"))
		maxDist = 1 - ::atof(cmdOpts.getOptStr("--identity"));

	if(cmdOpts.hasOpt("-s"))
		searchStrand = ::atoi(cmdOpts.getOptStr("-s"));
	if(cmdOpts.hasOpt("--strand"))
		searchStrand = ::atoi(cmdOpts.getOptStr("--strand"));

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* validate options */
	if(!(maxDist >= 0)) {
		cerr << "-i|--identity must between 0 and 1" << endl;
		return EXIT_FAILURE;
	}
	if(!(1 <= searchStrand && searchStrand <= 3)) {
		cerr << "-s|--strand must be 1, 2 or 3" << endl;
		return EXIT_FAILURE;
	}

	/* set filenames */
	msaFn = dbName + MSA_FILE_SUFFIX;
	csfmFn = dbName + CSFM_FILE_SUFFIX;
	hmmFn = dbName + HMM_FILE_SUFFIX;
	ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;

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

	seqIn.open(seqFn.c_str());
	if(!seqIn.is_open()) {
		cerr << "Unable to open seq file '" << seqFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	seqI.reset(&seqIn, AlphabetFactory::nuclAbc, seqFmt);

	/* open outputs */
	if(!outFn.empty()) {
		of.open(outFn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to '" << outFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}
	ostream& out = of.is_open() ? of : cout;

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
	size_t nNodes = ptu.numNodes();
	size_t nLeaves = ptu.numLeaves();
	const vector<PTUnrooted::PTUNodePtr>& id2node = ptu.getNodes();

	out << ANNEAL_HEADER << endl;
	const int K = hmm.getProfileSize();

	while(seqI.hasNext()) {
		PrimarySeq fwdRead = seqI.nextSeq();
		PrimarySeq revRead = fwdRead.revcom();
		string strand;
		double minCost = inf;
		int csStart = 0; /* 1-based csStart */
		int csEnd = 0; /* 1-based csEnd */
		string aln;
		BandedHMMP7::VScore seqVscore; /* construct an empty VScore */
		BandedHMMP7::VTrace seqVtrace; /* construct an empty VTrace */

		if(searchStrand & 01) {
			seqVscore.reset(fwdRead.length());
			hmm.calcViterbiScores(fwdRead, seqVscore); /* use original Viterbi algorithm */
			hmm.buildViterbiTrace(seqVscore, seqVtrace);
			if(seqVtrace.minScore < minCost) {
				strand = "+";
				csStart = hmm.getCSLoc(seqVtrace.alnStart);
				csEnd = hmm.getCSLoc(seqVtrace.alnEnd);
				aln = hmm.buildGlobalAlign(fwdRead, seqVscore, seqVtrace);
				minCost = seqVtrace.minScore;
			}
		}

		if(searchStrand & 02) {
			seqVscore.reset(revRead.length());
			hmm.calcViterbiScores(revRead, seqVscore); /* use original Viterbi algorithm */
			hmm.buildViterbiTrace(seqVscore, seqVtrace);
			if(seqVtrace.minScore < minCost) {
				strand = "-";
				csStart = hmm.getCSLoc(seqVtrace.alnStart);
				csEnd = hmm.getCSLoc(seqVtrace.alnEnd);
				aln = hmm.buildGlobalAlign(revRead, seqVscore, seqVtrace);
				minCost = seqVtrace.minScore;
			}
		}
		assert(aln.length() == csLen);
		assert(1 <= csStart && csStart <= csEnd && csEnd <= csLen);

		size_t hitsNodes = 0;
		size_t hitsLeaves = 0;

		for(vector<PTUnrooted::PTUNodePtr>::const_iterator node = id2node.begin(); node != id2node.end(); ++node) {
			double pDist = SeqUtils::pDist(aln, (*node)->getSeq(), csStart - 1, csEnd - 1); /* consider degenerated seq */
			if(pDist <= maxDist) {
				hitsNodes++;
				if((*node)->isLeaf())
					hitsLeaves++;
			}
		}

		/* output */
		out << fwdRead.getId() << "\t" << fwdRead.getDesc() << "\t" << fwdRead.getSeq() << "\t" <<
				strand << "\t" << csStart << "\t" << csEnd << "\t" << aln.substr(csStart - 1, csEnd - csStart + 1) << "\t" <<
				nNodes << "\t" << nLeaves << "\t" << hitsNodes << "\t" << hitsLeaves << "\t" <<
				static_cast<double>(hitsNodes) / nNodes << "\t" << static_cast<double>(hitsLeaves) / nLeaves << endl;
	}

	return 0;
}
