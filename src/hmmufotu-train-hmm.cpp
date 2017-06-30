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
 * Print the usage information of this program
 */
void printUsage(const string& progName) {
	cerr << "Train a Banded-HMM model used for HmmUFOtu program" << endl
		 << "Usage:    " << progName << "  <MSA-FILE> [options]" << endl
		 << "MSA-FILE  FILE                   : a multiple-alignment sequence file or pre-build MSA DB FILE" << endl
		 << "Options:    -o FILE              : write output to FILE instead of stdout" << endl
		 << "            -f|--symfrac DOUBLE  : conservation threshold for considering a site as a Match state in HMM [" << DEFAULT_SYMFRAC << "]" << endl
		 << "            -dm FILE             : use customized trained Dirichlet Model in FILE instead of the build-in file " << endl
		 << "            -v  FLAG             : enable verbose information" << endl
		 << "            -h|--help            : print this message and exit" << endl;
}

int main(int argc, char *argv[]) {
	ifstream in, dmin;
	ofstream of;
	double symfrac = DEFAULT_SYMFRAC;
	string infn;
	string outfn;
	string fmt;
	string dmfn;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.hasOpt("-h") || cmdOpts.hasOpt("--help")) {
		printUsage(argv[0]);
		return EXIT_SUCCESS;
	}
	if(cmdOpts.numMainOpts() != 1) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return EXIT_FAILURE;
	}

	infn = cmdOpts.getMainOpt(0);
	in.open(infn.c_str());
	if(!in.is_open()) {
		cerr << "Unable to open " << infn << " : " << ::strerror(errno)<< endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-o")) {
		outfn = cmdOpts.getOpt("-o");
		of.open(outfn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to " << outfn << endl;
			return EXIT_FAILURE;
		}
	}

	if(cmdOpts.hasOpt("-f"))
		symfrac = atof(cmdOpts.getOptStr("-f"));
	if(cmdOpts.hasOpt("-symfrac"))
		symfrac = atof(cmdOpts.getOptStr("-symfrac"));
	if(!(symfrac >= 0 && symfrac <= 1)) {
		cerr << "-f|--symfrac must between 0 and 1" << endl;
		return EXIT_FAILURE;
	}

	dmfn = PKG_DATADIR + string("/") + DEFAULT_DM_FILE;
	if(!ifstream(dmfn).good())
		dmfn = SRC_DATADIR + string("/") + DEFAULT_DM_FILE;

	if(cmdOpts.hasOpt("-dm"))
		dmfn = cmdOpts.getOpt("-dm");
	dmin.open(dmfn.c_str());
	if(!dmin.is_open()) {
		cerr << "Unable to open " << dmfn << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* guess input format */
	if(StringUtils::endsWith(infn, ".fasta") || StringUtils::endsWith(infn, ".fas")
		|| StringUtils::endsWith(infn, ".fa") || StringUtils::endsWith(infn, ".fna"))
		fmt = "fasta";
	else if(StringUtils::endsWith(infn, ".msa"))
		fmt = "msa";
	else {
		cerr << "Unrecognized MSA file format" << endl;
		return EXIT_FAILURE;
	}

	/* Load in BandedHmmPrior for the HMM training */
	BandedHMMP7Prior hmmPrior;
	dmin >> hmmPrior;
	if(dmin.bad()) {
		cerr << "Failed to read in the HMM Prior file " << endl;
		return EXIT_FAILURE;
	}

	/* Load data */
	MSA msa;
	if(fmt == "msa") { /* binary file provided */
		ifstream in(infn.c_str());
		msa.load(in);
		if(!in.good()) {
			cerr << "Unable to load MSA database" << endl;
			return EXIT_FAILURE;
		}
	}
	else {
		long nLoad = msa.loadMSAFile(ALPHABET, infn, fmt);
		if(nLoad <= 0) {
			cerr << "Unable to load MSA file" << endl;
			return EXIT_FAILURE;
		}
	}

	infoLog << "MSA loaded, found " << msa.getNumSeq() << " X " << msa.getCSLen() << " aligned sequences" << endl;

	ostream& out = outfn.empty() ? cout : of;

	BandedHMMP7 hmm; /* construct an empty profile */
	hmm.build(msa, symfrac, hmmPrior);
	infoLog << "Banded HMM profile trained" << endl;

	out << hmm;
	infoLog << "Banded HMM profile written" << endl;
}

