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
 * hmmufotu-sum.cpp
 *  summarize OTUTable from hmmufotu assignment files
 *  Created on: Apr 7, 2017
 *      Author: Qi Zheng
 *      Version: v1.1
 *
 */

#include <iostream>
#include <fstream>
#include <strstream>
#include <cctype>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <limits>
#include <map>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp> /* for boost string join */
#include <boost/lexical_cast.hpp>
#include "HmmUFOtu.h"
#include "HmmUFOtu_main.h"

using namespace std;
using namespace EGriceLab;
using namespace Eigen;

/* default values */
static const string ALIGN_FORMAT = "fasta";
static const double DEFAULT_EFFN = 2;
static const int DEFAULT_MIN_NREAD = 0;
static const int DEFAULT_MIN_NSAMPLE = 0;
static const double DEFAULT_MIN_Q = 0;
static const double DEFAULT_MIN_ALN_IDENTITY = 0;
static const double DEFAULT_MIN_HMM_IDENTITY = 0;
typedef boost::unordered_map<PTUnrooted::PTUNodePtr, OTUObserved> OTUMap;
typedef boost::unordered_map<string, vector<string>> OTU2ReadMap;

/**
 * Print introduction of this program
 */
void printIntro(void) {
	cerr << "Generate phylogeny-based OTUTable with taxonomy information, and optional consensus based OTU representatives and tree" << endl;
}

/**
 * Print the usage information
 */
void printUsage(const string& progName) {
	cerr << "Usage:    " << progName << "  <HmmUFOtu-DB> <(INFILE [INFILE2 ...]> <-o OTU-OUT> [options]" << endl
		 << "INFILE          FILE           : assignment file(s) from hmmufotu, by default the filenames will be used as sample names" << endl
		 << "Options:    -o  FILE           : OTU summary output, required" << endl
		 << "            -r  FILE           : output the read IDs for each OTU" << endl
		 << "            -l  FILE           : an optional sample list, with 1st field sample-name and 2nd field input filename, only input files in list will be used" << endl
		 << "            -c  FILE           : OTU Consensus Sequence (CS) alignment of each OTU" << endl
		 << "            -t  FILE           : OTU tree output" << endl
		 << "            -q  DBL            : minimum qTaxon score (negative log10 posterior error rate) required [" << DEFAULT_MIN_Q << "]" << endl
		 << "            --aln-iden DBL     : minimum alignment identity required for assignment result [" << DEFAULT_MIN_ALN_IDENTITY << "]" << endl
		 << "            --hmm-iden DBL     : minimum profile-HMM identity required for assignment result [" << DEFAULT_MIN_HMM_IDENTITY << "]" << endl
		 << "            -e|--effN  DBL     : effective number of sequences (pseudo-count) for inferring CS of OTUs with Dirichelet Density models, set 0 to disable [" << DEFAULT_EFFN << "]" << endl
		 << "            -n  INT            : minimum number of observed reads required to define an OTU across all samples, 0 for no filtering [" << DEFAULT_MIN_NREAD << "]" << endl
		 << "            -s  INT            : minimum number of observed samples required to define an OTU, 0 for no filtering [" << DEFAULT_MIN_NSAMPLE << "]" << endl
		 << "            --no-gap  FLAG     : if -c is set, this will output the non-gapped OTU sequences instead of aligned CS alignment" << endl
		 << "            -v  FLAG           : enable verbose information, you may set multiple -v for more details" << endl
		 << "            --version          : show program version and exit" << endl
		 << "            -h|--help          : print this message and exit" << endl;
}

int main(int argc, char* argv[]) {
	/* variable declarations */
	string dbName, msaFn, hmmFn, ptuFn;
	vector<string> inFiles;
	map<string, string> sampleFn2Name;
	string listFn;
	string otuFn, readFn, csFn, treeFn;
	ifstream msaIn, hmmIn, ptuIn;
	ofstream otuOut, readOut, treeOut, csOut;
	SeqIO csO;
	OTU2ReadMap otu2Read;

	double effN = DEFAULT_EFFN;
	double minQ = DEFAULT_MIN_Q;
	double minAlnIden = DEFAULT_MIN_ALN_IDENTITY;
	double minHmmIden = DEFAULT_MIN_HMM_IDENTITY;
	int minRead = DEFAULT_MIN_NREAD;
	int minSample = DEFAULT_MIN_NSAMPLE;
	bool noGap = false;

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
		sampleFn2Name[fn] = fn; /* use filename as samplename by default */
	}

	if(cmdOpts.hasOpt("-o"))
		otuFn = cmdOpts.getOpt("-o");
	else {
		cerr << "-o must be specified" << endl;
		return EXIT_FAILURE;
	}
	if(cmdOpts.hasOpt("-r"))
		readFn = cmdOpts.getOpt("-r");
	if(cmdOpts.hasOpt("-c"))
		csFn = cmdOpts.getOpt("-c");
	if(cmdOpts.hasOpt("-t"))
		treeFn = cmdOpts.getOpt("-t");

	if(cmdOpts.hasOpt("-l"))
		listFn = cmdOpts.getOpt("-l");

	if(cmdOpts.hasOpt("-e"))
		effN = ::atof(cmdOpts.getOptStr("-e"));
	if(cmdOpts.hasOpt("--effN"))
		effN = ::atof(cmdOpts.getOptStr("--effN"));

	if(cmdOpts.hasOpt("-q"))
		minQ = ::atof(cmdOpts.getOptStr("-q"));
	if(cmdOpts.hasOpt("--aln-iden"))
		minAlnIden = ::atof(cmdOpts.getOptStr("--aln-iden"));
	if(cmdOpts.hasOpt("--hmm-iden"))
		minHmmIden = ::atof(cmdOpts.getOptStr("--hmm-iden"));

	if(cmdOpts.hasOpt("-n"))
		minRead = ::atoi(cmdOpts.getOptStr("-n"));
	if(cmdOpts.hasOpt("-s"))
		minSample = ::atoi(cmdOpts.getOptStr("-s"));

	if(cmdOpts.hasOpt("--no-gap"))
		noGap = true;

	if(cmdOpts.hasOpt("-v"))
		INCREASE_LEVEL(cmdOpts.getOpt("-v").length());

	/* validate options */
	if(!(effN >= 0)) {
		cerr << "-e|--effN must be non-negative" << endl;
		return EXIT_FAILURE;
	}
	if(!(minRead >= 0)) {
		cerr << "-n must be non-negative integer" << endl;
		return EXIT_FAILURE;
	}
	if(!(minSample >= 0)) {
		cerr << "-s must be non-negative integer" << endl;
		return EXIT_FAILURE;
	}

	/* set filenames */
	msaFn = dbName + MSA_FILE_SUFFIX;
	hmmFn = dbName + HMM_FILE_SUFFIX;
	ptuFn = dbName + PHYLOTREE_FILE_SUFFIX;

	/* open inputs */
	if(!listFn.empty()) {
		ifstream listIn(listFn.c_str());
		int nRead = 0;
		if(!listIn.is_open()) {
			cerr << "Unable to open sample list '" << listFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		infoLog << "Read in sample names from " << listFn << endl;
		inFiles.clear(); /* clear inFiles */
//		sampleFn2Name.clear(); /* clear sample names */
		string line;
		while(std::getline(listIn, line)) {
			if(line[0] == '#')
				continue;
			vector<string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			if(fields.size() >= 2) {
				string name = fields[0];
				string fn = fields[1];
				if(sampleFn2Name.count(fn) > 0) { /* this is an input file */
					inFiles.push_back(fn);
					sampleFn2Name[fn] = name; /* update the sample name */
					nRead++;
				}
			}
		}
		listIn.close();
		infoLog << nRead << " user-provided sample names read" << endl;
	}

	msaIn.open(msaFn.c_str(), ios_base::in | ios_base::binary);
	if(!msaIn) {
		cerr << "Unable to open MSA data '" << msaFn << "': " << ::strerror(errno) << endl;
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

	/* open outputs */
	otuOut.open(otuFn.c_str());
	if(!otuOut.is_open()) {
		cerr << "Unable to write to '" << otuFn << "': " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(!readFn.empty()) {
		readOut.open(readFn.c_str());
		if(!readOut.is_open()) {
			cerr << "Unable to write to '" << readFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
	}

	if(!csFn.empty()) {
		csOut.open(csFn.c_str());
		if(!csOut.is_open()) {
			cerr << "Unable to write to '" << csFn << "': " << ::strerror(errno) << endl;
			return EXIT_FAILURE;
		}
		csO.reset(&csOut, AlphabetFactory::nuclAbc, ALIGN_FORMAT);
	}

	if(!treeFn.empty()) {
		treeOut.open(treeFn.c_str());
		if(!treeOut.is_open()) {
			cerr << "Unable to write to '" << treeFn << "': " << ::strerror(errno) << endl;
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
	ptu.setRoot(0);

	const DegenAlphabet* abc = msa.getAbc();
	const int S = inFiles.size();
	const int L = ptu.numAlignSites();
	const size_t N = ptu.numNodes();

	/* process input files */
	OTUMap otuData;
	vector<string> sampleNames;
	for(int s = 0; s < inFiles.size(); ++s) {
		string infn = inFiles[s];
		string sample = sampleFn2Name[infn];
		infoLog << "Processing sample " << sampleFn2Name[infn] << " ..." << endl;
		ifstream in(infn.c_str());
		if(!in.is_open()) {
			cerr << "Unable to open '" << infn << "': " << ::strerror(errno) << " sampled ignored" << endl;
			continue;
		}
		TSVScanner tsvIn(in, true);
		sampleNames.push_back(sample);
		while(tsvIn.hasNext()) {
			const TSVRecord& record = tsvIn.nextRecord();
			string rid = record.getFieldByName("id");
			int csStart = ::atoi(record.getFieldByName("CS_start").c_str());
			int csEnd = ::atoi(record.getFieldByName("CS_end").c_str());
			const string& aln = record.getFieldByName("alignment");
			const long taxon_id = ::atol(record.getFieldByName("taxon_id").c_str());
			double qTaxon = ::atof(record.getFieldByName("Q_taxon").c_str());

			if(taxon_id >= 0 && qTaxon >= minQ
					&& EGriceLab::alignIdentity(abc, aln, csStart - 1, csEnd -1)
					&& EGriceLab::hmmIdentity(hmm, aln, csStart - 1, csEnd - 1)) { /* a valid assignment */
				const PTUnrooted::PTUNodePtr& node = ptu.getNode(taxon_id);
				string otuID = boost::lexical_cast<string>(node->getId());
				if(otuData.count(node) == 0) /* not initiated */
					otuData.insert(std::make_pair(node, OTUObserved(otuID, node->getTaxon(), L, S)));
				OTUObserved& otu = otuData.find(node)->second;
				otu.count(s)++;
				if(readOut.is_open())
					otu2Read[otuID].push_back(rid);
				for(int j = 0; j < L; ++j) {
					int8_t b = abc->encode(::toupper(aln[j]));
					if(b >= 0)
						otu.freq(b, j)++;
					else
						otu.gap(j)++;
				}
			}
		}
	}

	/* construct an OTU table and output alignment */
	infoLog << "Computing OTUTable" << endl;
	OTUTable otuTable(sampleNames);

	for(int i = 0; i < N; ++i) {
		const PTUnrooted::PTUNodePtr& node = ptu.getNode(i);
		if(otuData.count(node) == 0) // not an observed OTU
			continue;
		OTUObserved& otu = otuData.find(node)->second;
		if(otu.numReads() >= minRead && otu.numSamples() >= minSample) /* filter OTUs */
			otuTable.addOTU(otu);
		else {/* remove unnessesary otu2Read elements */
			if(readOut.is_open())
				otu2Read.erase(otu.id);
		}
	}

	/* write the OTU table */
	infoLog << "Writing OTUTable" << endl;
	otuOut << otuTable;

	/* write read list */
	if(readOut.is_open()) {
		infoLog <<"Wring read ID lists" << endl;
		for(size_t i = 0; i <otuTable.numOTUs(); ++i) {
			string otuID = otuTable.getOTU(i);
			readOut << otuID <<"\t" << boost::join(otu2Read[otuID], " ") << endl;
		}
	}

	/* write the CS seq */
	if(csOut.is_open()) {
		infoLog << "Writing OTU Consensus Sequences" << endl;
		for(size_t i = 0; i < otuTable.numOTUs(); ++i) {
			PTUnrooted::PTUNodePtr node = ptu.getNode(::atol(otuTable.getOTU(i).c_str()));
			OTUObserved& data = otuData.find(node)->second;
			int nRead = data.count.sum();
			int nSample = (data.count.array() > 0).count();

			DigitalSeq csSeq = ptu.inferPostCS(node, data.freq, data.gap, effN);
			string desc = "DBName="
					+ dbName + ";Taxonomy=\"" + otuTable.getTaxon(i)
					+ "\";AnnoDist=" + boost::lexical_cast<string>(node->getAnnoDist())
					+ ";ReadCount=" + boost::lexical_cast<string>(nRead)
					+ ";SampleHits=" + boost::lexical_cast<string>(nSample);
			PrimarySeq otuSeq(csSeq.getAbc(), csSeq.getName(), csSeq.toString(), desc);
			if(noGap)
				otuSeq.removeGaps();
			csO.writeSeq(otuSeq);
		}
	}

	/* write the tree */
	if(treeOut.is_open()) {
		infoLog << "Writing OTU tree" << endl;

		boost::unordered_set<PTUnrooted::PTUNodePtr> otuSet;
		for(size_t i = 0; i < otuTable.numOTUs(); ++i)
			otuSet.insert(ptu.getNode(boost::lexical_cast<long> (otuTable.getOTU(i))));

		ptu.exportTree(treeOut, otuSet, "newick");
	}
}
