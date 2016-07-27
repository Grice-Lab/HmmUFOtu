/*
 * hmmufotu-build-dm.cpp
 *
 *  Created on: Jul 15, 2016
 *      Author: zhengqi
 */

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include "HmmUFOtu.h"

using namespace std;
using namespace Eigen;
using namespace EGriceLab;
using namespace EGriceLab::Math;

static const int DEFAULT_QM = 1;
static const double DEFAULT_SYMFRAC = 0.5;

/**
 * Print the usage information of this program
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <MSA-INFILE> [options]" << endl
		 << "Options:    -o FILE    : write output to FILE instead of stdout" << endl
		 << "            -qM INT    : number of Dirichlet Mixture model components for match state emissions [" << DEFAULT_QM << "]" << endl
		 << "            -symfrac   : conservation threshold for an MSA site to be considered as a Match state [" << DEFAULT_SYMFRAC << "]" << endl;
}

static const int MAX_NUM_COMPO = 4;

int main(int argc, char* argv[]) {
	ifstream in;
	filebuf* fb; /* only be used when -o specified */
	ostream out(NULL);
	int qM = DEFAULT_QM;
	double symfrac = DEFAULT_SYMFRAC;
	string infn;
	string fmt;

	/* parse options */
	CommandOptions cmdOpts(argc, argv);
	if(cmdOpts.numMainOpts() != 1) {
		cerr << "Error:" << endl;
		printUsage(argv[0]);
		return -1;
	}
	infn = cmdOpts.getMainOpt(0);
	in.open(infn.c_str());
	if(!in.is_open()) {
		cerr << "Unable to open '" << cmdOpts.getMainOpt(0) << "'" << endl;
		return -1;
	}

	if(cmdOpts.hasOpt("-o")) {
		fb->open(cmdOpts.getOpt("-o").c_str(), ios::out);
		if(!fb->is_open()) {
			cerr << "Unable to write to '" << cmdOpts.getOpt("-o") << "'" << endl;
			return -1;
		}
		out.rdbuf(fb);
	}
	else
		out.rdbuf(cout.rdbuf());

	if(cmdOpts.hasOpt("-qM"))
		qM = atoi(cmdOpts.getOpt("-qM").c_str());
	if(!(qM > 0 && qM < MAX_NUM_COMPO)) {
		cerr << "-qM must between 1 and " << MAX_NUM_COMPO << endl;
		return -1;
	}

	if(cmdOpts.hasOpt("-symfrac"))
		symfrac = atof(cmdOpts.getOpt("-symfrac").c_str());
	if(!(symfrac >= 0 && symfrac <= 1)) {
		cerr << "-symfrac must between 0 and 1" << endl;
	}

	/* guess input format */
	if(StringUtils::endsWith(infn, ".fasta") || StringUtils::endsWith(infn, ".fas")
		|| StringUtils::endsWith(infn, ".fa") || StringUtils::endsWith(infn, ".fna"))
		fmt = "fasta";
	else if(StringUtils::endsWith(infn, ".msa"))
		fmt = "msa";
	else {
		cerr << "Unrecognized MSA file format" << endl;
		return -1;
	}

	/* Load data */
	MSA* msa;
	if(fmt == "msa") { /* binary file provided */
		ifstream in(infn.c_str());
		msa = MSA::load(in);
	}
	else
		msa = MSA::loadMSAFile("dna", infn, fmt); /* always read DNA MSA file */
	cerr << "MSA loaded" << endl;
	msa->prune(); /* prune MSA */
	cerr << "MSA pruned" << endl;

	const DegenAlphabet* abc = msa->getAbc();
	const int K = abc->getSize();
	assert (K == 4);
	/* construct DM models */
	DirichletModel* dmME;
	if(qM == 1)
		dmME = new DirichletDensity(K);
	else
		dmME = new DirichletMixture(K, qM); /* Match state emission, 4 categories, qM components */
	DirichletModel* dmIE = new DirichletDensity(K); /* Insert state emission */
	DirichletModel* dmMT = new DirichletDensity(3); /* Match state transition, 3 categories */
	DirichletModel* dmIT = new DirichletDensity(2); /* Insert state transition, 3 categories */
	DirichletModel* dmDT = new DirichletDensity(2); /* Delete state transition, 3 categories */

	cerr << "Dirichlet models initiated" << endl;

	/* Prepare training data */
	Matrix4Xd dataME(K, msa->getCSLen());
	Matrix4Xd dataIE(K, msa->getCSLen());
	Matrix3Xd dataMT = Matrix3Xd::Zero(3, msa->getCSLen()); /* M->M, M->I and M->D */
	Matrix2Xd dataIT = Matrix2Xd::Zero(2, msa->getCSLen()); /* I->M, I->I and I->D */
	Matrix2Xd dataDT = Matrix2Xd::Zero(2, msa->getCSLen()); /* D->M, D->I and D->D */

	int cME = 0;
	int cIE = 0;
	int cMT = 0;
	int cIT = 0;
	int cDT = 0;

	for(int j = 0; j < msa->getCSLen(); ++j) {
		if(msa->symWFrac(j) >= symfrac) /* match state emission */
			dataME.col(cME++) = msa->symWFreq(j);
		else
			dataIE.col(cIE++) = msa->symWFreq(j);
	}
	dataME.conservativeResize(K, cME);
	dataIE.conservativeResize(K, cIE);
	cerr << "Emission training data prepared" << endl;

	const unsigned L = msa->getCSLen();
	for(int j = 0; j < L - 1; ++j) {
		bool matchFlag = msa->symWFrac(j) >= symfrac;
		for(int i = 0; i < msa->getNumSeq(); ++i) {
			double w = msa->getSeqWeight(i);
			bool resFlag = msa->encodeAt(i, j) >= 0;
			if(!matchFlag && !resFlag) /* ignore phentome positions */
				continue;
			/* search to next non-phentome position */
			bool matchFlagN = false;
			bool resFlagN = false;
			int k = j + 1;
			while(!matchFlagN && !resFlagN && k < L) {
				matchFlagN = msa->symWFrac(k) >= symfrac;
				resFlagN = msa->encodeAt(i, k) >= 0;
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
		}
		if((dataMT.col(cMT).array() != 0).any())
			cMT++;
		if((dataIT.col(cIT).array() != 0).any())
			cIT++;
		if((dataDT.col(cDT).array() != 0).any())
			cDT++;
	}
	cerr << "cMT:" << cMT << " cIT:" << cIT << " cDT:" << cDT << endl;
	dataMT.conservativeResize(3, cMT);
	dataIT.conservativeResize(2, cIT);
	dataDT.conservativeResize(2, cDT);
	cerr << "Transition training data prepared" << endl;

	/* train DM models */
	double costME = dmME->trainML(dataME);
	cerr << "Match emission model trained" << endl;
	double costIE = dmIE->trainML(dataIE);
	cerr << "Insert emission model trained" << endl;
	double costMT = dmMT->trainML(dataMT);
	cerr << "Match transition model trained" << endl;
	double costIT = dmIT->trainML(dataIT);
	cerr << "Insert transition model trained" << endl;
	double costDT = dmDT->trainML(dataDT);
	cerr << "Delete transition model trained" << endl;

	/* output */
	out << "Match emission:" << endl << *dmME << endl;
	out << "Insert emission:" << endl << *dmIE << endl;
	out << "Match transition:" << endl << *dmMT << endl;
	out << "Insert transition:" << endl << *dmIT << endl;
	out << "Delete transition:" << endl << *dmDT << endl;
}
