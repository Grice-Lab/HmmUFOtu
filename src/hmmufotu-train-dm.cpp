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
 * hmmufotu-build-dm.cpp
 * train a customized Dirichlet Model using a MSA file or database
 *  Created on: Jul 15, 2016
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <math.h> /* use C99 header */
#include <boost/iostreams/filtering_stream.hpp> /* basic boost streams */
#include <boost/iostreams/device/file.hpp> /* file sink and source */
#include <boost/iostreams/filter/zlib.hpp> /* for zlib support */
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp> /* for bzip2 support */
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_hmm.h"
#include "EGMath.h"

using namespace std;
using namespace Eigen;
using namespace EGriceLab;
using namespace EGriceLab::Math;
using namespace EGriceLab::HmmUFOtu;

static const int DEFAULT_QM = 5;
static const double DEFAULT_SYMFRAC = 0.5;
static const int MAX_NUM_COMPO = 10;
static const double DEFAULT_PRI_RATE = 0.05;
static const int MAX_ITER = 0;
static const int DEFAULT_NSEED = 1;
static const string ALPHABET = "dna";

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Train an HmmUFOtu prior model using Dirichlet Density/Mixture models with customized data" << endl;
}

/**
 * Print the usage information of this program
 */
void printUsage(const string& progName) {
	string ZLIB_SUPPORT;
	#ifdef HAVE_LIBZ
	ZLIB_SUPPORT = ", support .gz or .bz2 compressed file";
	#endif

	cerr << "Usage:    " << progName << "  <MSA-FILE> [options]" << endl
		 << "MSA-FILE  FILE             : a multiple-alignment sequence file or pre-build MSA DB FILE" << ZLIB_SUPPORT << endl
		 << "Options:    -o FILE        : write output to FILE instead of stdout" << endl
		 << "            --fmt  STR     : MSA format, supported format: 'fasta', 'msa'" << endl
		 << "            -qM INT[>=2]   : number of Dirichlet Mixture model components for match state emissions [" << DEFAULT_QM << "]" << endl
		 << "            -symfrac       : conservation threshold for an MSA site to be considered as a Match state [" << DEFAULT_SYMFRAC << "]" << endl
		 << "            --max-it INT   : maximum iteration allowed in gradient descent training, 0 for no limit [" << MAX_ITER << "]" << endl
		 << "            --pri-rate DBL : adjust the sequence weights so the prior information is roughly this ratio in training [" << DEFAULT_PRI_RATE << "]" << endl
		 << "            -s|--seed INT  : random seed used in Dirichlet Mixture model training (-qM > 1) for debug purpose" << endl
		 << "            -n  INT        : number of different random seeds in Dirichlet Mixture model training [" << DEFAULT_NSEED << "]" << endl
		 << "            -v  FLAG       : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version      : show program version and exit" << endl
		 << "            -h|--help      : print this help and exit" << endl;
}

int main(int argc, char* argv[]) {
	boost::iostreams::filtering_istream in;
	ofstream of;
	int qM = DEFAULT_QM;
	double symfrac = DEFAULT_SYMFRAC;
	double priRate = DEFAULT_PRI_RATE;
	int maxIter = MAX_ITER;
	string inFn;
	string outFn;
	string fmt;
	unsigned seed = time(NULL); // using time as default seed
	int nSeed = DEFAULT_NSEED;

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

	if(cmdOpts.hasOpt("-qM"))
		qM = ::atoi(cmdOpts.getOpt("-qM").c_str());
	if(!(qM > 1 && qM <= MAX_NUM_COMPO)) {
		cerr << "-qM must between 2 and " << MAX_NUM_COMPO << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-symfrac"))
		symfrac = ::atof(cmdOpts.getOpt("-symfrac").c_str());
	if(!(symfrac >= 0 && symfrac <= 1)) {
		cerr << "-symfrac must between 0 and 1" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("--pri-rate"))
		priRate = ::atof(cmdOpts.getOpt("--pri-rate").c_str());
	if(!( priRate > 0 && priRate <= 1 )) {
		cerr << "--rate must be in (0, 1]" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("--max-it"))
		maxIter = ::atoi(cmdOpts.getOpt("--max-it").c_str());
	if(maxIter < 0) {
		cerr << "--max-it must be a non-negative integer" << endl;
		return EXIT_FAILURE;
	}

	if(cmdOpts.hasOpt("-s"))
		seed = ::atoi(cmdOpts.getOpt("-s").c_str());
	if(cmdOpts.hasOpt("--seed"))
		seed = ::atoi(cmdOpts.getOpt("--seed").c_str());

	if(cmdOpts.hasOpt("-n"))
		nSeed = ::atoi(cmdOpts.getOpt("-n").c_str());

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* guess input format */
	if(fmt.empty()) {
		if(StringUtils::endsWith(inFn, ".msa"))
			fmt = "msa";
		else {
			string inPre = inFn;
			StringUtils::removeEnd(inPre, GZIP_FILE_SUFFIX);
			StringUtils::removeEnd(inPre, BZIP2_FILE_SUFFIX);
			fmt = SeqUtils::guessSeqFileFormat(inPre);
		}
	}
	if(!(fmt == "fasta" || fmt == "msa")) {
		cerr << "Unsupported sequence format '" << fmt << "'" << endl;
		return EXIT_FAILURE;
	}

	/* set random seed */
	srand(seed);

	/* open input */
#ifdef HAVE_LIBZ
	if(StringUtils::endsWith(inFn, GZIP_FILE_SUFFIX))
		in.push(boost::iostreams::gzip_decompressor());
	else if(StringUtils::endsWith(inFn, BZIP2_FILE_SUFFIX))
		in.push(boost::iostreams::bzip2_decompressor());
	else { }
#endif
	/* open source */
	boost::iostreams::file_source inSrc(inFn);
	if(!inSrc.is_open()) {
		cerr << "Unable to open seq file '" << inFn << "' " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	in.push(inSrc);

	/* open output */
	if(!outFn.empty()) {
		of.open(outFn.c_str());
		if(!of.is_open()) {
			cerr << "Unable to write to '" << outFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}
	ostream& out = of.is_open() ? of : cout;

	/* Load data */
	MSA msa;
	if(fmt == "msa") { /* binary file provided */
		if(loadProgInfo(in).bad())
			return EXIT_FAILURE;
		msa.load(in);
	}
	else {
		msa.loadMSA(ALPHABET, in, fmt);
		msa.setName(inFn);
	}

	if(!in.bad()) /* load sequence format */
		infoLog << "MSA loaded" << endl;
	else {
		cerr << "Unable to load MSA seq from '" << inFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	if(!msa.pruned()) {
		msa.prune(); /* prune MSA if necessary*/
		infoLog << "MSA pruned" << endl;
	}
	infoLog << "MSA database created for " << msa.getNumSeq() << " X " << msa.getCSLen() << " aligned sequences" << endl;

	double effN = 1 / priRate;
	msa.sclaleWeight(effN / msa.getNumSeq());
	infoLog << "MSA total weight scaled as: " << effN << endl;

	const int K = msa.getAbc()->getSize();
	assert (K == 4);
	/* construct an HMM prior */
	BandedHMMP7Prior pri;
	/* set the # of parameters */
	pri.setDims(K, qM);
	pri.setMaxIter(maxIter);

	infoLog << "Dirichlet prior model initiated" << endl;

	const unsigned L = msa.getCSLen();
	const unsigned N = msa.getNumSeq();
	/* Prepare the training data */
	Matrix4Xd dataME(K, L);
	Matrix4Xd dataIE(K, L);
	Matrix3Xd dataMT = Matrix3Xd::Zero(3, L); /* M->M, M->I and M->D */
	Matrix2Xd dataIT = Matrix2Xd::Zero(2, L); /* I->M and I->I */
	Matrix2Xd dataDT = Matrix2Xd::Zero(2, L); /* D->M and D->D */

	int cME = 0;
	int cIE = 0;
	int cMT = 0;
	int cIT = 0;
	int cDT = 0;

	for(int j = 0; j < L; ++j) {
		if(msa.symWFrac(j) >= symfrac) /* match state emission */
			dataME.col(cME++) = msa.symWFreq(j);
		else
			dataIE.col(cIE++) = msa.symWFreq(j);
	}
	dataME.conservativeResize(K, cME);
	dataIE.conservativeResize(K, cIE);
//	cerr << "Emission training data prepared" << endl;
//	cerr << "cME:" << cME << endl;
//	cerr << "cIE:" << cIE << endl;

	for(int j = 0; j < L - 1; ++j) {
		bool matchFlag = msa.symWFrac(j) >= symfrac;
		for(int i = 0; i < N; ++i) {
			double w = msa.getSeqWeight(i);
			bool resFlag = msa.encodeAt(i, j) >= 0;
			if(!matchFlag && !resFlag) /* ignore phantom positions */
				continue;
			/* search to next non-phantom position */
			bool matchFlagN = false;
			bool resFlagN = false;
			int k = j + 1;
			while(!matchFlagN && !resFlagN && k < L) {
				matchFlagN = msa.symWFrac(k) >= symfrac;
				resFlagN = msa.encodeAt(i, k) >= 0;
				k++;
			}
			if(k >= L) /* no more position found */
				continue;
			/* Match state transition */
			if(matchFlag && resFlag) {
				if(matchFlagN && resFlagN) /* M->M */
					dataMT(0, cMT) += w;
				else if(!matchFlagN && resFlagN) /* M->I */
					dataMT(1, cMT) += w;
				else if(matchFlagN && !resFlagN) /* M->D */
					dataMT(2, cMT) += w;
				else { }
			}
			/* Insert state transition */
			else if(!matchFlag && resFlag) {
				if(matchFlagN && resFlagN) /* I->M */
					dataIT(0, cIT) += w;
				else if(!matchFlagN && resFlagN) /* I->I */
					dataIT(1, cIT) += w;
				else { }
			}
			/* Delete state transition */
			else if(matchFlag && !resFlag) {
				if(matchFlagN && resFlagN) /* D->M */
					dataDT(0, cDT) += w;
				else if(matchFlagN && !resFlagN) /* D->D */
					dataDT(1, cDT) += w;
				else { } // ignore other cases
			}
			else { }
		} /* end each seq */
		if((dataMT.col(cMT).array() != 0).any())
			cMT++;
		if((dataIT.col(cIT).array() != 0).any())
			cIT++;
		if((dataDT.col(cDT).array() != 0).any())
			cDT++;
	} /* end each position */

	dataMT.conservativeResize(3, cMT);
	dataIT.conservativeResize(2, cIT);
	dataDT.conservativeResize(2, cDT);
	infoLog << "Transition training data prepared" << endl;

	/* train DM models */
	/* iteratively train ME */
	double costME = inf;
	int bestIdx = 0;

	/* make a copy of the original model */
	DirichletMixture model(pri.dmME);
	infoLog << "Training Match Emission model" <<endl;
	for(int i = 1; i <= nSeed; ++i) {
		double cost = model.trainML(dataME);
		cerr << "  seed " << i << " trained, cost: " << cost << endl;
		if(cost < costME) { // a better model found
			pri.dmME = model; // copy back
			bestIdx = i;
			costME = cost;
		}
	}
	if(!::isnan(costME))
		infoLog << "Best Match Emission model found at seed " << bestIdx << endl;
	else {
		cerr << "Unable to train Match Emission model" << endl;
		return EXIT_FAILURE;
	}

	double costIE = pri.dmIE.trainML(dataIE);
	infoLog << "Insert Emission model trained" << endl;

	double costMT = pri.dmMT.trainML(dataMT);
	infoLog << "Match Transition model trained" << endl;

	double costIT = pri.dmIT.trainML(dataIT);
	infoLog << "Insert Transition model trained" << endl;

	double costDT = pri.dmDT.trainML(dataDT);
	infoLog << "Delete Transition model trained" << endl;

	/* output */
	out << pri;
}
