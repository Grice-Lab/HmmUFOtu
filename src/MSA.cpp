/*
 * MSA.cpp
 *
 *  Created on: Jul 23, 2015
 *      Author: zhengqi
 */

#include <fstream>
#include <map>
#include <cstdlib>
#include <cctype>
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

double MSA::identityAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // check range once
		throw out_of_range("CS pos is out of range");
	return resCount.col(j).maxCoeff() / static_cast<double>(numSeq);
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
	vector<unsigned> pruningSites;
	for(unsigned j = 0; j < csLen; ++j) {
		unsigned numGap = 0;
		for(unsigned i = 0; i < numSeq; ++i)
			if(abc->isGap(concatMSA[i * csLen + j]))
				numGap++;
		if(numGap == numSeq)
			pruningSites.push_back(j);
	}
	if(pruningSites.empty()) /* nothing to do */
		return *this;
	/* pruning concatMSA backward */
	for(unsigned i = numSeq; i > 0; i--)
		for(vector<unsigned>::const_reverse_iterator rit = pruningSites.rbegin(); rit != pruningSites.rend(); ++rit)
			concatMSA.erase((i-1) * csLen + *rit, 1);
	/* pruning the known CS, if exist */
	if(!CS.empty())
		for(vector<unsigned>::const_reverse_iterator rit = pruningSites.rbegin(); rit != pruningSites.rend(); ++rit)
			CS.erase(*rit, 1);
	/* rebuild the counts */
	updateRawCounts();
	updateSeqWeight();
	updateWeightedCounts();
	isPruned = true;
	return *this;
}

MSA* MSA::loadFastaFile(const string& alphabet, const string& filename) {
	MSA* msa = new MSA(alphabet);
	msa->setName(basename(filename));
	SeqIO seqIn(filename, alphabet, "fasta");
	while(seqIn.hasNext()) {
		const PrimarySeq& seq = seqIn.nextSeq();
		//cerr << seq.getId() << " " << seq.getSeq() << endl;
		/* check new seq */
		if(msa->csLen != 0 && seq.length() != msa->csLen)
			throw ios_base::failure("Not a valid fasta alignment file where not all sequences have the same length!");
		msa->csLen = seq.length(); /* update csLen */
		msa->numSeq++;
		msa->totalNumGap += seq.numGap();
		//ids.push_back(seq.getId());
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
	delete[] resCountBuf;
	delete[] gapCountBuf;
	delete[] seqWeightBuf;
	delete[] resWCountBuf;
	delete[] gapWCountBuf;
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
		if(resCount.col(j).maxCoeff() >= gapCount(j)) { /* consensus is not a gap */
			MatrixXi::Index maxRow;
			resCount.col(j).maxCoeff(&maxRow);
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
	for(unsigned i = 0; i < numSeq; ++i)
		for(unsigned j = 0; j < csLen; ++j) {
			char c = ::toupper(residualAt(i, j));
			if(abc->isSymbol(c))
				resCount(abc->encode(c), j)++;
			else if(abc->isGap(c))
				gapCount(j)++;
			else { } // do nothing
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
		seqWeight(i) = w;
	}
	/* normalize in place */
	seqWeight /= numSeq;
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
				gapWCount(j) += seqWeight(j);
			else { } // do nothing
		}
}

bool MSA::save(ostream& f) {
	/* get aux length */
	string::size_type nAlphabet = alphabet.length();
	string::size_type nCS = CS.length();
	string::size_type nName = name.length();
	/* write basic info */
	f.write((const char*) &nAlphabet, sizeof(string::size_type));
	f.write(alphabet.c_str(), nAlphabet + 1);
	f.write((const char*) &nName, sizeof(string::size_type));
	f.write(name.c_str(), nName + 1);
//	cerr << "alphabet written" << endl;
	f.write((const char*) &numSeq, sizeof(unsigned));
	f.write((const char*) &csLen, sizeof(unsigned));
	f.write((const char*) &totalNumGap, sizeof(unsigned long));
	f.write((const char*) &nCS, sizeof(string::size_type));
	f.write(CS.c_str(), nCS + 1); /* write the null terminal */
	f.write((const char*) &isPruned, sizeof(bool));
//	cerr << "CS written" << endl;
	/* write concatMSA */
	f.write(concatMSA.c_str(), concatMSA.length() + 1);
//	cerr << "concatMSA written" << endl;
	/* write raw counts */
	f.write((const char*) resCountBuf, sizeof(int) * abc->getSize() * csLen);
	f.write((const char*) gapCountBuf, sizeof(int) * csLen);
	/* write seq weights */
	f.write((const char*) seqWeightBuf, sizeof(double) * numSeq);
	/* write weighted counts */
	f.write((const char*) resWCountBuf, sizeof(double) * abc->getSize() * csLen);
	f.write((const char*) gapWCountBuf, sizeof(double) * csLen);
	return f.good();
}


MSA* MSA::load(istream& f) {
	/* construct a MSA object */
	MSA* msa = new MSA(); // use default alphabet first
	/* read basic info */
	char* buf = NULL;
	string::size_type nAlphabet, nCS, nName;
	f.read((char*) &nAlphabet, sizeof(string::size_type));
	buf = new char[nAlphabet + 1];
	f.read(buf, nAlphabet + 1); /* read the null terminal */
	msa->alphabet.assign(buf, nAlphabet); // override the original alphabet
	delete[] buf;
	msa->abc = SeqCommons::getAlphabetByName(msa->alphabet);
	f.read((char*) &nName, sizeof(string::size_type));
	buf = new char[nName + 1];
	f.read(buf, nName + 1); /* read the null terminal */
	msa->name.assign(buf, nName);
	delete[] buf;

	f.read((char*) &msa->numSeq, sizeof(unsigned));
	f.read((char*) &msa->csLen, sizeof(unsigned));
	f.read((char*) &msa->totalNumGap, sizeof(unsigned long));
	f.read((char*) &nCS, sizeof(string::size_type));
	buf = new char[nCS + 1];
	f.read(buf, nCS + 1);
	msa->CS.assign(buf, nCS); /* read the null terminal */
	delete[] buf;
	f.read((char*) &msa->isPruned, sizeof(bool));

	/* read concatMSA */
	buf = new char[msa->numSeq * msa->csLen + 1];
	f.read(buf, msa->numSeq * msa->csLen + 1);
	msa->concatMSA.assign(buf, msa->numSeq * msa->csLen);
	delete[] buf;
	/* Read raw counts */
	msa->resCountBuf = new int[msa->abc->getSize() * msa->csLen];
	f.read((char*) msa->resCountBuf, sizeof(int) * msa->abc->getSize() * msa->csLen);
	msa->gapCountBuf = new int[msa->csLen];
	f.read((char*) msa->gapCountBuf, sizeof(int) * msa->csLen);
	/* Read seq weights */
	msa->seqWeightBuf = new double[msa->numSeq];
	f.read((char*) msa->seqWeightBuf, sizeof(double) * msa->numSeq);
	/* Read weighted counts */
	msa->resWCountBuf = new double[msa->abc->getSize() * msa->csLen];
	f.read((char*) msa->resWCountBuf, sizeof(double) * msa->abc->getSize() * msa->csLen);
	msa->gapWCountBuf = new double[msa->csLen];
	f.read((char*) msa->gapWCountBuf, sizeof(double) * msa->csLen);

	/* replacement constructor the count matrices */
	new (&msa->resCount) Map<MatrixXi>(msa->resCountBuf, msa->abc->getSize(), msa->csLen);
	new (&msa->gapCount) Map<VectorXi>(msa->gapCountBuf, msa->csLen);
	new (&msa->seqWeight) Map<VectorXd>(msa->seqWeightBuf, msa->numSeq);
	new (&msa->resWCount) Map<MatrixXd>(msa->resWCountBuf, msa->abc->getSize(), msa->csLen);
	new (&msa->gapWCount) Map<VectorXd>(msa->gapWCountBuf, msa->csLen);

	return msa;
}

} /* namespace EGriceLab */
