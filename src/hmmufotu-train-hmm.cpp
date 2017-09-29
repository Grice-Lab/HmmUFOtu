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
 * hmmufotu-train.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: zhengqi
 */


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <ctime>
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_hmm.h"
#include "EGMath.h"

#ifndef PKG_DATADIR
#define PKG_DATADIR "."
#endif

#ifndef SRC_DATADIR
#define SRC_DATADIR "."
#endif

using namespace std;
using namespace Eigen;
using namespace EGriceLab;
using namespace EGriceLab::Math;

static const double DEFAULT_SYMFRAC = 0.5;
static const string DEFAULT_DM_FILE = "gg_97_otus.dm";
static const string ALPHABET = "dna";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Train a Banded-HMM model with customized data" << endl;
}

/**
 * Print the usage information of this program
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <MSA-FILE> [options]" << endl
		 << "MSA-FILE  FILE                   : a multiple-alignment sequence file or pre-build MSA DB FILE" << endl
		 << "Options:    -o FILE              : write output to FILE instead of stdout" << endl
		 << "            --fmt  STR           : MSA format, supported format: 'fasta', 'msa'" << endl
		 << "            -f|--symfrac DOUBLE  : conservation threshold for considering a site as a Match state in HMM [" << DEFAULT_SYMFRAC << "]" << endl
		 << "            -dm FILE             : use customized trained Dirichlet Model in FILE instead of the build-in file " << endl
		 << "            -v  FLAG             : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version            : show program version and exit" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char *argv[]) {
	ifstream in, dmIn;
	ofstream of;
	double symfrac = DEFAULT_SYMFRAC;
	string inFn;
	string outFn;
	string fmt;
	string dmFn;

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

	if(cmdOpts.numMainOpts() != 1) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	inFn = cmdOpts.getMainOpt(0);


	if(cmdOpts.hasOpt("-o"))
		outFn = cmdOpts.getOpt("-o");

	if(cmdOpts.hasOpt("--fmt"))
		fmt = cmdOpts.getOpt("--fmt");

	if(cmdOpts.hasOpt("-f"))
		symfrac = atof(cmdOpts.getOptStr("-f"));
	if(cmdOpts.hasOpt("-symfrac"))
		symfrac = atof(cmdOpts.getOptStr("-symfrac"));
	if(!(symfrac >= 0 && symfrac <= 1)) {
		cerr << "-f|--symfrac must between 0 and 1" << endl;
		return EXIT_FAILURE;
	}

	dmFn = PKG_DATADIR + string("/") + DEFAULT_DM_FILE;
	if(!ifstream(dmFn.c_str()).good())
		dmFn = SRC_DATADIR + string("/") + DEFAULT_DM_FILE;

	if(cmdOpts.hasOpt("-dm"))
		dmFn = cmdOpts.getOpt("-dm");
	dmIn.open(dmFn.c_str());
	if(!dmIn.is_open()) {
		cerr << "Unable to open " << dmFn << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* guess input format */
	if(fmt.empty()) {
		if(StringUtils::endsWith(inFn, ".fasta") || StringUtils::endsWith(inFn, ".fas")
		|| StringUtils::endsWith(inFn, ".fa") || StringUtils::endsWith(inFn, ".fna"))
			fmt = "fasta";
		else if(StringUtils::endsWith(inFn, ".msa"))
			fmt = "msa";
		else {
			cerr << "Unrecognized MSA file format" << endl;
			return EXIT_FAILURE;
		}
	}
	if(!(fmt == "fasta" || fmt == "msa")) {
		cerr << "Unsupported sequence format '" << fmt << "'" << endl;
		return EXIT_FAILURE;
	}

	/* open input and output */
	if(fmt != "msa")
		in.open(inFn.c_str());
	else
		in.open(inFn.c_str(), ios_base::binary);
	if(!in.is_open()) {
		cerr << "Unable to open " << inFn << " : " << ::strerror(errno)<< endl;
		return EXIT_FAILURE;
	}

	if(!outFn.empty()) {
		of.open(outFn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to " << outFn << endl;
			return EXIT_FAILURE;
		}
	}
	/* Load in BandedHmmPrior for the HMM training */
	BandedHMMP7Prior hmmPrior;
	dmIn >> hmmPrior;
	if(dmIn.bad()) {
		cerr << "Failed to read in the HMM Prior file " << endl;
		return EXIT_FAILURE;
	}

	/* Load data */
	MSA msa;
	long nLoad;
	if(fmt == "msa") /* binary file provided */
		msa.load(in);
	else
		nLoad = msa.loadMSA(ALPHABET, in, fmt);
	if(!in.bad() && nLoad >= 0) /* load sequence format */
		infoLog << "MSA loaded" << endl;
	else {
		cerr << "Unable to load MSA seq from '" << inFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(!msa.pruned()) {
		msa.prune(); /* prune MSA if necessary*/
		infoLog << "MSA pruned" << endl;
	}
	infoLog << "MSA database created for " << msa.getNumSeq() << " X " << msa.getCSLen() << " aligned sequences" << endl;

	ostream& out = outFn.empty() ? cout : of;

	BandedHMMP7 hmm; /* construct an empty profile */
	hmm.build(msa, symfrac, hmmPrior);
	infoLog << "Banded HMM profile trained" << endl;

	out << hmm;
	infoLog << "Banded HMM profile written" << endl;
}

