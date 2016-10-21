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

char MSA::CSResidualAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // check range once
		throw out_of_range("CS pos is out of range");

	return CS.empty() ? '\0' /* not calculated yet */ : CS[j];
}

char MSA::CSBaseAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // check range once
		throw out_of_range("CS pos is out of range");

	MatrixXd::Index max;
	VectorXd freq = resWCount.col(j);
	freq.maxCoeff(&max);
	return abc->decode(max);
}

double MSA::identityAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // check range once
		throw out_of_range("CS pos is out of range");
	return resCount.col(j).maxCoeff() / static_cast<double>(numSeq);
}

double MSA::wIdentityAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // check range once
		throw out_of_range("CS pos is out of range");
	return resWCount.col(j).maxCoeff() / numSeq;
}

double MSA::gapFrac(unsigned j) const {
	return gapCount(j) / static_cast<double>(numSeq);
}

double MSA::symFrac(unsigned j) const {
	double numRes = resCount.col(j).sum();
	double numGap = gapCount(j);
	return numRes / (numRes + numGap);
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

	/* rebuild the counts */
	updateRawCounts();
	updateSeqWeight();
	updateWeightedCounts();
	isPruned = true;
	return *this;
}

MSA* MSA::loadFastaFile(const string& alphabet, const string& filename) {
	MSA* msa = new MSA(alphabet);
	msa->setName(StringUtils::basename(filename));
	SeqIO seqIn(filename, alphabet, "fasta");
	while(seqIn.hasNext()) {
		const PrimarySeq& seq = seqIn.nextSeq();
		//cerr << seq.getId() << " " << seq.getSeq() << endl;
		/* check new seq */
		if(msa->csLen != 0 && seq.length() != msa->csLen)
			throw ios_base::failure("Not a valid fasta alignment file where not all sequences have the same length!");
		msa->csLen = seq.length(); /* update csLen */
		msa->numSeq++;
		//ids.push_back(seq.getId());
		msa->seqNames.push_back(seq.getId());
		msa->concatMSA.append(seq.getSeq());
	}
	assert(msa->concatMSA.length() == msa->numSeq * msa->csLen);
	msa->updateRawCounts();
	msa->updateSeqWeight();
	msa->updateWeightedCounts();
	msa->calculateCS();
	return msa;
}

void MSA::clear() {
	delete[] startIdx;
	delete[] endIdx;
	delete[] lenIdx;
	delete[] resCountBuf;
	delete[] gapCountBuf;
	delete[] seqWeightBuf;
	delete[] resWCountBuf;
	delete[] gapWCountBuf;
	startIdx = NULL;
	endIdx = NULL;
	lenIdx = NULL;
	resCountBuf = NULL;
	gapCountBuf = NULL;
	seqWeightBuf = NULL;
	resWCountBuf = NULL;
	gapWCountBuf = NULL;
}

void MSA::resetRawCount() {
	/* Initiate raw buffers and associated them to the Eigen maps */
	if(resCountBuf == NULL) {
		resCountBuf = new int[abc->getSize() * csLen];
		new (&resCount) Map<MatrixXi>(resCountBuf, abc->getSize(), csLen); /* replacement constructor */
	}
	if(gapCountBuf == NULL) {
		gapCountBuf = new int[csLen];
		new (&gapCount) Map<VectorXi>(gapCountBuf, csLen); /* replacement constructor */
	}
	if(startIdx == NULL)
		startIdx = new int[numSeq](); /* zero initiation */
	if(endIdx == NULL)
		endIdx = new int[numSeq]; /* zero initiation */
	if(lenIdx == NULL)
		lenIdx = new int[numSeq]; /* zero initiation */
	/* reset count to zero */
	resCount.setZero();
	gapCount.setZero();
}

void MSA::resetSeqWeight() {
	/* Initiate raw buffers and associated them to the Eigen maps */
	if(seqWeightBuf == NULL) {
		seqWeightBuf = new double[numSeq];
		new (&seqWeight) Map<VectorXd>(seqWeightBuf, numSeq); /* replacement constructor */
	}
	/* reset count to zero */
	seqWeight.setZero();
}

void MSA::resetWeightedCount() {
	/* Initiate raw buffers and associated them to the Eigen maps */
	if(resWCountBuf == NULL) {
		resWCountBuf = new double[abc->getSize() * csLen];
		new (&resWCount) Map<MatrixXd>(resWCountBuf, abc->getSize(), csLen); /* replacement constructor */
	}
	if(gapWCountBuf == NULL) {
		gapWCountBuf = new double[csLen];
		new (&gapWCount) Map<VectorXd>(gapWCountBuf, csLen); /* replacement constructor */
	}
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
	/* write auxiliary index */
	out.write((const char*) startIdx, sizeof(int) * numSeq);
	out.write((const char*) endIdx, sizeof(int) * numSeq);
	out.write((const char*) lenIdx, sizeof(int) * numSeq);

	/* write raw counts */
	out.write((const char*) resCountBuf, sizeof(int) * abc->getSize() * csLen);
	out.write((const char*) gapCountBuf, sizeof(int) * csLen);
	/* write seq weights */
	out.write((const char*) seqWeightBuf, sizeof(double) * numSeq);
	/* write weighted counts */
	out.write((const char*) resWCountBuf, sizeof(double) * abc->getSize() * csLen);
	out.write((const char*) gapWCountBuf, sizeof(double) * csLen);
	return out;
}

MSA* MSA::load(istream& in) {
	/* construct a MSA object */
	MSA* msa = new MSA(); // use default alphabet first
	/* Read program info */
	string pname, pver;
	readProgName(in, pname);
	if(pname != progName) {
		cerr << "Not a MSA object file" << endl;
		in.setstate(ios_base::failbit);
		return msa;
	}
	readProgVersion(in, pver);
	if(cmpVersion(progVersion, pver) < 0) {
		cerr << "You are trying using an older version " << (progName + progVersion) <<
				" to read a newer MSA object file that was build by " << (pname + pver) << endl;
		in.setstate(ios_base::failbit);
		return msa;
	}
	/* read basic info */
	char* buf = NULL;
	string::size_type nAlphabet, nCS, nName;
	in.read((char*) &nAlphabet, sizeof(string::size_type));
	buf = new char[nAlphabet + 1];
	in.read(buf, nAlphabet + 1); /* read the null terminal */
	msa->alphabet.assign(buf, nAlphabet); // override the original alphabet
	delete[] buf;
	msa->abc = SeqCommons::getAlphabetByName(msa->alphabet);
	in.read((char*) &nName, sizeof(string::size_type));
	buf = new char[nName + 1];
	in.read(buf, nName + 1); /* read the null terminal */
	msa->name.assign(buf, nName);
	delete[] buf;

	in.read((char*) &msa->numSeq, sizeof(unsigned));
	in.read((char*) &msa->csLen, sizeof(unsigned));
	in.read((char*) &nCS, sizeof(string::size_type));
	buf = new char[nCS + 1];
	in.read(buf, nCS + 1);
	msa->CS.assign(buf, nCS); /* read the null terminal */
	delete[] buf;
	in.read((char*) &msa->isPruned, sizeof(bool));
	/* read seqNames */
	msa->seqNames.resize(msa->numSeq); /* set all names to empty */
	for(vector<string>::iterator it = msa->seqNames.begin(); it != msa->seqNames.end();) {
		char c = in.get();
		if(c != '\0')
			it->push_back(c);
		else
			it++;
	}
	/* read concatMSA */
	buf = new char[msa->numSeq * msa->csLen + 1];
	in.read(buf, msa->numSeq * msa->csLen + 1);
	msa->concatMSA.assign(buf, msa->numSeq * msa->csLen);
	delete[] buf;
	/* Read auxiliary index */
	msa->startIdx = new int[msa->numSeq];
	in.read((char*) msa->startIdx, sizeof(int) * msa->numSeq);
	msa->endIdx = new int[msa->numSeq];
	in.read((char*) msa->endIdx, sizeof(int) * msa->numSeq);
	msa->lenIdx = new int[msa->numSeq];
	in.read((char*) msa->lenIdx, sizeof(int) * msa->numSeq);
	/* Read raw counts */
	msa->resCountBuf = new int[msa->abc->getSize() * msa->csLen];
	in.read((char*) msa->resCountBuf, sizeof(int) * msa->abc->getSize() * msa->csLen);
	msa->gapCountBuf = new int[msa->csLen];
	in.read((char*) msa->gapCountBuf, sizeof(int) * msa->csLen);
	/* Read seq weights */
	msa->seqWeightBuf = new double[msa->numSeq];
	in.read((char*) msa->seqWeightBuf, sizeof(double) * msa->numSeq);
	/* Read weighted counts */
	msa->resWCountBuf = new double[msa->abc->getSize() * msa->csLen];
	in.read((char*) msa->resWCountBuf, sizeof(double) * msa->abc->getSize() * msa->csLen);
	msa->gapWCountBuf = new double[msa->csLen];
	in.read((char*) msa->gapWCountBuf, sizeof(double) * msa->csLen);

	/* replacement constructor the count matrices */
	new (&msa->resCount) Map<MatrixXi>(msa->resCountBuf, msa->abc->getSize(), msa->csLen);
	new (&msa->gapCount) Map<VectorXi>(msa->gapCountBuf, msa->csLen);
	new (&msa->seqWeight) Map<VectorXd>(msa->seqWeightBuf, msa->numSeq);
	new (&msa->resWCount) Map<MatrixXd>(msa->resWCountBuf, msa->abc->getSize(), msa->csLen);
	new (&msa->gapWCount) Map<VectorXd>(msa->gapWCountBuf, msa->csLen);

	return msa;
}

} /* namespace EGriceLab */
