/*
 * MSA.cpp
 *
 *  Created on: Jul 23, 2015
 *      Author: zhengqi
 */

#include <fstream>
#include <cstdlib>
#include <cctype>
#include <set>
#include <algorithm>
#include "HmmUFOtuConst.h"
#include "MSA.h"
#include "Stats.h"
#include "LinearAlgebraBasic.h"
#include "SeqIO.h"
#include "StringUtils.h"

namespace EGriceLab {
using namespace std;
using namespace Math;

/* static field definition */
const double DEFAULT_CONSENSUS_FRAC = 0.5;

char MSA::CSResidualAt(unsigned j) const {
	return CS.empty() ? '\0' /* not calculated yet */ : CS.at(j);
}

char MSA::CSBaseAt(unsigned j) const {
	MatrixXd::Index max;
	VectorXd freq = resWCount.col(j);
	freq.maxCoeff(&max);
	return abc->decode(max);
}

double MSA::identityAt(unsigned j) const {
	return resCount.col(j).maxCoeff() / static_cast<double>(numSeq);
}

double MSA::wIdentityAt(unsigned j) const {
	return resWCount.col(j).maxCoeff() / numSeq;
}

double MSA::gapFrac(unsigned j) const {
	return gapCount(j) / static_cast<double>(numSeq);
}

double MSA::gapWFrac(unsigned j) const {
	double numRes = resWCount.col(j).sum();
	double numGap = gapWCount(j);
	return numGap / (numRes + numGap);
}

double MSA::symFrac(unsigned j) const {
	return resCount.col(j).sum() / static_cast<double>(numSeq);
}

double MSA::symWFrac(unsigned j) const {
	double numRes = resWCount.col(j).sum();
	double numGap = gapWCount(j);
	return numRes / (numRes + numGap);
}

MSA& MSA::prune() {
	if(isPruned)
		return *this;
	if(numSeq == 0)
		return *this;
	set<unsigned> pruningSites;
	for(unsigned j = 0; j < csLen; ++j)
		if(resCount.col(j).sum() == 0) // no residual at this site
			pruningSites.insert(j);

	if(pruningSites.empty()) /* nothing to do */
		return *this;

	/* construct the pruned concatMSA */
	string prunedMSA;
	prunedMSA.reserve(concatMSA.length() - numSeq * pruningSites.size());
	/* copy the concatMSA to prunedMSA, ignore pruned sites */
	for(string::size_type i = 0; i != concatMSA.length(); ++i)
		if(pruningSites.find(i % csLen) == pruningSites.end()) // this is not a pruning site
			prunedMSA.push_back(concatMSA[i]);
	/* swap the storage */
	concatMSA.swap(prunedMSA);

	/* pruning the known CS, if exist */
	if(!CS.empty()) {
		string prunedCS;
		prunedCS.reserve(CS.length() - pruningSites.size());
		/* copy the CS to prunedCS, ignore pruned sites */
		for(string::size_type j = 0; j != CS.length(); ++j)
			if(pruningSites.find(j) == pruningSites.end()) // this is not a pruning site
				prunedCS.push_back(CS[j]);
		/* swap the storage */
		CS.swap(prunedCS);
	}

	/* update index */
	csLen -= pruningSites.size();

	/* destroy old counts */
	clear();
	resetRawCount();
	resetSeqWeight();
	resetWeightedCount();

	/* rebuild the counts */
	updateRawCounts();
	updateSeqWeight();
	updateWeightedCounts();
	isPruned = true;

	return *this;
}

long MSA::loadFastaFile(const string& alphabet, const string& filename) {
	setName(StringUtils::basename(filename));
	SeqIO seqIn(filename, alphabet, "fasta");
	while(seqIn.hasNext()) {
		const PrimarySeq& seq = seqIn.nextSeq();
		//cerr << seq.getId() << " " << seq.getSeq() << endl;
		/* check new seq */
		if(csLen != 0 && seq.length() != csLen) {
			cerr << "Invalid fasta alignment file! Not all sequences have the same length!";
			return -1;
		}
		csLen = seq.length(); /* update csLen */
		numSeq++;
		//ids.push_back(seq.getId());
		seqNames.push_back(seq.getId());
		concatMSA.append(seq.getSeq());
	}
	assert(concatMSA.length() == numSeq * csLen);
	updateRawCounts();
	updateSeqWeight();
	updateWeightedCounts();
	calculateCS();

	return numSeq;
}

void MSA::clear() {
//	resCount.setZero();
//	gapCount.setZero();
//	seqWeight.setZero();
//	resWCount.setZero();
//	gapWCount.setZero();
}

void MSA::resetRawCount() {
	/* Initiate counts */
	if(resCount.size() != abc->getSize() * csLen)
		resCount.resize(abc->getSize(), csLen);
	if(gapCount.rows() != csLen)
		gapCount.resize(csLen);
	if(startIdx.size() != numSeq)
		startIdx.resize(numSeq);
	if(endIdx.size() != numSeq)
		endIdx.resize(numSeq);
	if(lenIdx.size() != numSeq)
		lenIdx.resize(numSeq);

	/* reset count to zero */
	resCount.setZero();
	gapCount.setZero();
	std::fill(startIdx.begin(), startIdx.end(), 0);
	std::fill(endIdx.begin(), endIdx.end(), 0);
	std::fill(lenIdx.begin(), lenIdx.end(), 0);
}

void MSA::resetSeqWeight() {
	if(seqWeight.rows() != numSeq)
		seqWeight.resize(numSeq);
	/* reset count to zero */
	seqWeight.setZero();
}

void MSA::resetWeightedCount() {
	if(resWCount.size() != abc->getSize() * csLen)
		resWCount.resize(abc->getSize(), csLen);
	if(gapWCount.rows() != csLen)
		gapWCount.resize(csLen);
	/* reset count to zero */
	resWCount.setZero();
	gapWCount.setZero();
}

void MSA::calculateCS() {
	if(CS.length() == csLen) /* already calculated, ignore */
		return;
	CS.clear();
	for(unsigned j = 0; j < csLen; ++j) {
		char csRes;
		if(resWCount.col(j).maxCoeff() >= gapWCount(j)) { /* consensus is not a gap */
			MatrixXi::Index maxRow;
			resWCount.col(j).maxCoeff(&maxRow);
			csRes = abc->decode(maxRow);
		}
		else
		    csRes = abc->getGap()[0]; /* use the first gap character */
		CS.push_back(csRes);
	}
}

void MSA::updateRawCounts() {
	/* reset old data */
	resetRawCount();
	/* calculate raw count */
	for(int i = 0; i < numSeq; ++i) {
		int start = -1;
		int end = -1;
		int len = 0;
		for(int j = 0; j < csLen; ++j) {
			char c = ::toupper(residualAt(i, j));
			if(abc->isSymbol(c)) {
				if(start == -1)
					start = j;
				if(j > end)
					end = j;
				len++;
				resCount(abc->encode(c), j)++;
			}
			else if(abc->isGap(c))
				gapCount(j)++;
			else { } // do nothing
		}
		startIdx[i] = start;
		endIdx[i] = end;
		lenIdx[i] = len;
	}
}

void MSA::updateSeqWeight() {
	/* reset old data */
	resetSeqWeight();
	/* Get a pssw weight matrix */
	MatrixXi pssw(abc->getSize(), csLen);
	for(unsigned j = 0; j != csLen; ++j)
		pssw.col(j) = (resCount.col(j).array() != 0).count() * resCount.col(j);

	/* get seq weights by summing over all CS pos */
	for(unsigned i = 0; i != numSeq; ++i) {
		double w = 0;
		for(unsigned j = 0; j != csLen; ++j) {
			int8_t b = encodeAt(i, j);
			if(b >= 0) /* is a valid symbol */
				w += 1.0 / pssw(b, j);
		}
		if(seqLength(i) > 0)
			w /= seqLength(i); /* first normalize weight by non-gap seqLength */
		seqWeight(i) = w;
	}
	/* bring seqWeight to nseq */
	seqWeight *= numSeq / seqWeight.sum();
}

void MSA::updateWeightedCounts() {
	/* reset old data */
	resetWeightedCount();
	/* calculate weighted count */
	for(unsigned i = 0; i < numSeq; ++i)
		for(unsigned j = 0; j < csLen; ++j) {
			char c = ::toupper(residualAt(i, j));
			if(abc->isSymbol(c))
				resWCount(abc->encode(c), j) += seqWeight(i);
			else if(abc->isGap(c))
				gapWCount(j) += seqWeight(i);
			else { } // do nothing
		}
}

ostream& MSA::save(ostream& out) {
	/* save program info */
	writeProgName(out, progName);
	writeProgVersion(out, progVersion);
	/* get aux length */
	string::size_type nAlphabet = alphabet.length();
	string::size_type nCS = CS.length();
	string::size_type nName = name.length();
	/* write basic info */
	out.write((const char*) &nAlphabet, sizeof(string::size_type));
	out.write(alphabet.c_str(), nAlphabet + 1);
	out.write((const char*) &nName, sizeof(string::size_type));
	out.write(name.c_str(), nName + 1);
//	cerr << "alphabet written" << endl;
	out.write((const char*) &numSeq, sizeof(unsigned));
	out.write((const char*) &csLen, sizeof(unsigned));
	out.write((const char*) &nCS, sizeof(string::size_type));
	out.write(CS.c_str(), nCS + 1); /* write the null terminal */
	//	cerr << "CS written" << endl;
	out.write((const char*) &isPruned, sizeof(bool));
	/* write seqNames */
	for(vector<string>::const_iterator it = seqNames.begin(); it != seqNames.end(); ++it)
//		f.write((const char*) it->length(), sizeof(string::size_type)); /* write seq length first */
		out.write(it->c_str(), it->length() + 1); /* write the null terminal */
	/* write concatMSA */
	out.write(concatMSA.c_str(), concatMSA.length() + 1);
//	cerr << "concatMSA written" << endl;

	int* bufi = NULL; /* integer output buffer */
	double* bufd = NULL; /* double output buffer */

	/* write auxiliary index */
	bufi = new int[numSeq];
	std::copy(startIdx.begin(), startIdx.end(), bufi);
	out.write((const char*) bufi, sizeof(int) * numSeq);

	std::copy(endIdx.begin(), endIdx.end(), bufi);
	out.write((const char*) bufi, sizeof(int) * numSeq);

	std::copy(lenIdx.begin(), lenIdx.end(), bufi);
	out.write((const char*) bufi, sizeof(int) * numSeq);
	delete[] bufi;

	/* write raw counts */
	bufi = new int[resCount.size()];
	Map<MatrixXi> resCountMap(bufi, resCount.rows(), resCount.cols());
	resCountMap = resCount; /* copy data */
	out.write((const char*) bufi, sizeof(int) * resCount.size());
	delete[] bufi;

	bufi = new int[gapCount.rows()];
	Map<VectorXi> gapCountMap(bufi, gapCount.rows());
	gapCountMap = gapCount; /* copy data */
	out.write((const char*) bufi, sizeof(int) * gapCount.rows());
	delete[] bufi;

	/* write seq weights */
	bufd = new double[seqWeight.rows()];
	Map<VectorXd> seqWeightMap(bufd, seqWeight.rows());
	seqWeightMap = seqWeight; /* copy data */
	out.write((const char*) bufd, sizeof(double) * seqWeight.rows());
	delete[] bufd;

	/* write weighted counts */
	bufd = new double[resWCount.size()];
	Map<MatrixXd> resWCountMap(bufd, resWCount.rows(), resWCount.cols());
	resWCountMap = resWCount; /* copy data */
	out.write((const char*) bufd, sizeof(double) * resWCount.size());
	delete[] bufd;

	bufd = new double[gapWCount.rows()];
	Map<VectorXd> gapWCountMap(bufd, gapWCount.rows());
	gapWCountMap = gapWCount; /* copy data */
	out.write((const char*) bufd, sizeof(double) * gapWCount.rows());
	delete[] bufd;

	return out;
}

istream& MSA::load(istream& in) {
	/* Read program info */
	string pname, pver;
	readProgName(in, pname);
	if(pname != progName) {
		cerr << "Not a MSA object file" << endl;
		in.setstate(ios_base::failbit);
		return in;
	}
	readProgVersion(in, pver);
	if(cmpVersion(progVersion, pver) < 0) {
		cerr << "You are trying using an older version " << (progName + progVersion) <<
				" to read a newer MSA data file that was build by " << (pname + pver) << endl;
		in.setstate(ios_base::failbit);
		return in;
	}
	/* read basic info */
	char* buf = NULL; /* character buf */
	int* bufi = NULL; /* integer buf */
	double* bufd = NULL; /* double buf */
	string::size_type nAlphabet, nCS, nName;
	in.read((char*) &nAlphabet, sizeof(string::size_type));
	buf = new char[nAlphabet + 1];
	in.read(buf, nAlphabet + 1); /* read the null terminal */
	alphabet.assign(buf, nAlphabet); // override the original alphabet
	delete[] buf;
	abc = AlphabetFactory::getAlphabetByName(alphabet);
	in.read((char*) &nName, sizeof(string::size_type));
	buf = new char[nName + 1];
	in.read(buf, nName + 1); /* read the null terminal */
	name.assign(buf, nName);
	delete[] buf;

	in.read((char*) &numSeq, sizeof(unsigned));
	in.read((char*) &csLen, sizeof(unsigned));
	in.read((char*) &nCS, sizeof(string::size_type));
	buf = new char[nCS + 1];
	in.read(buf, nCS + 1);
	CS.assign(buf, nCS); /* read the null terminal */
	delete[] buf;
	in.read((char*) &isPruned, sizeof(bool));

	/* read seqNames */
	seqNames.resize(numSeq); /* set all names to empty */
	for(vector<string>::iterator it = seqNames.begin(); it != seqNames.end();) {
		char c = in.get();
		if(c != '\0')
			it->push_back(c);
		else
			it++;
	}

	/* read concatMSA */
	buf = new char[numSeq * csLen + 1];
	in.read(buf, numSeq * csLen + 1);
	concatMSA.assign(buf, numSeq * csLen);
	delete[] buf;

	/* initiate all maticies and indices */
	resetRawCount();
	resetSeqWeight();
	resetWeightedCount();

	/* Read auxiliary index */
	bufi = new int[numSeq];
	in.read((char*) bufi, sizeof(int) * numSeq);
	std::copy(bufi, bufi + numSeq, startIdx.begin()); /* copy data */

	in.read((char*) bufi, sizeof(int) * numSeq);
	std::copy(bufi, bufi + numSeq, endIdx.begin()); /* copy data */

	in.read((char*) bufi, sizeof(int) * numSeq);
	std::copy(bufi, bufi + numSeq, lenIdx.begin()); /* copy data */
	delete[] bufi;

	/* Read raw counts */
	bufi = new int[resCount.size()];
	in.read((char*) bufi, sizeof(int) * resCount.size());
	resCount = Map<MatrixXi>(bufi, resCount.rows(), resCount.cols()); /* copy by assign */
	delete[] bufi;

	bufi = new int[gapCount.rows()];
	in.read((char*) bufi, sizeof(int) * gapCount.rows());
	gapCount = Map<VectorXi>(bufi, gapCount.rows()); /* copy by assign */
	delete[] bufi;

	/* Read seq weights */
	bufd = new double[seqWeight.rows()];
	in.read((char*) bufd, sizeof(double) * seqWeight.rows()); /* copy by assign */
	seqWeight = Map<VectorXd>(bufd, seqWeight.rows());
	delete[] bufd;

	/* Read weighted counts */
	bufd = new double[resWCount.size()];
	in.read((char*) bufd, sizeof(double) * resWCount.size());
	resWCount = Map<MatrixXd>(bufd, resWCount.rows(), resWCount.cols()); /* copy by assign */
	delete[] bufd;

	bufd = new double[gapWCount.rows()];
	in.read((char*) bufd, sizeof(double) * gapWCount.rows());
	gapWCount = Map<VectorXd>(bufd, gapWCount.rows()); /* copy by assign */
	delete[] bufd;

	return in;
}

} /* namespace EGriceLab */
