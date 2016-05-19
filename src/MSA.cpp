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
#include "SeqIO.h"
#include "StringUtils.h"

namespace EGriceLab {
using namespace std;

char MSA::CSResidualAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // check range once
		throw out_of_range("CS pos is out of range");

	return CS.empty() ? '\0' /* not calculated yet */ : CS[j];
}

double MSA::identityAt(unsigned j) const {
	if(!(j >= 0 && j < csLen)) // check range once
		throw out_of_range("CS pos is out of range");
	return max(resCount[j], abc->getSize()) / static_cast<double>(numSeq);
}

double MSA::percentGapAt(unsigned j) const {
	return gapCount[j] / static_cast<double>(numSeq);
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
	updateCounts();
	isPruned = true;
	return *this;
}

MSA* MSA::loadFastaFile(const string& alphabet, const string& filename) {
	MSA* msa = new MSA(alphabet);
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
	msa->updateCounts();

	msa->calculateCS();
	return msa;
}

void MSA::clear() {
	if(resCount != NULL) {
		for(unsigned j = 0; j != csLen; ++j)
			delete[] resCount[j]; /* delete inner arrays */
		delete[] resCount;
	}
	delete[] gapCount;
}

void MSA::calculateCS() {
	if(CS.length() == csLen) /* already calculated, ignore */
		return;
	CS.clear();
	for(unsigned j = 0; j < csLen; ++j) {
		char csRes;
		if(max(resCount[j], abc->getSize()) >= gapCount[j]) /* consensus is not a gap */
			csRes = abc->decode(which_max(resCount[j], abc->getSize()));
		else
		    csRes =  abc->getGap()[0]; /* use the first gap character */
		CS.push_back(csRes);
	}
}

void MSA::updateCounts() {
	/* delete old data */
	clear();
	/* Build the count matrices */
	resCount = new unsigned*[csLen];
	for(unsigned j = 0; j != csLen; ++j)
		resCount[j] = new unsigned[abc->getSize()](); /* zero initiation for inner arrays */
	gapCount = new unsigned[csLen](); /* use zero initiation */

	for(unsigned i = 0; i < numSeq; ++i)
		for(unsigned j = 0; j < csLen; ++j) {
			char c = ::toupper(residualAt(i, j));
			if(abc->isSymbol(c))
				resCount[j][abc->encode(c)]++;
			else if(abc->isGap(c))
				gapCount[j]++;
			else { } // do nothing
		}
}

bool MSA::save(ostream& f) {
	/* write basic info */
	string::size_type nAlphabet = alphabet.length();
	string::size_type nCS = CS.length();
	f.write((const char*) &nAlphabet, sizeof(string::size_type));
	f.write(alphabet.c_str(), nAlphabet + 1);
//	cerr << "alphabet written" << endl;
	f.write((const char*) &numSeq, sizeof(unsigned));
	f.write((const char*) &csLen, sizeof(unsigned));
	f.write((const char*) &totalNumGap, sizeof(unsigned long));
	f.write((const char*) &nCS, sizeof(string::size_type));
	f.write(CS.c_str(), nCS + 1); /* write the null terminal */
//	cerr << "CS written" << endl;
	/* write concatMSA */
	f.write(concatMSA.c_str(), concatMSA.length() + 1);
//	cerr << "concatMSA written" << endl;
	/* write count matrices */
	for(unsigned j = 0; j < csLen; ++j)
		f.write((const char*) resCount[j], sizeof(unsigned) * abc->getSize());
	f.write((const char*) gapCount, sizeof(unsigned) * csLen);
	return f.good();
}


MSA* MSA::load(istream& f) {
	/* construct a MSA object */
	MSA* msa = new MSA(); // use default alphabet first
	/* read basic info */
	char* buf = NULL;
	string::size_type nAlphabet, nCS;
	f.read((char*) &nAlphabet, sizeof(string::size_type));
	buf = new char[nAlphabet + 1];
	f.read(buf, nAlphabet + 1); /* read the null terminal */
	msa->alphabet.assign(buf, nAlphabet); // override the original alphabet
	delete[] buf;
	msa->abc = SeqCommons::getAlphabetByName(msa->alphabet);
	f.read((char*) &msa->numSeq, sizeof(unsigned));
	f.read((char*) &msa->csLen, sizeof(unsigned));
	f.read((char*) &msa->totalNumGap, sizeof(unsigned long));
	f.read((char*) &nCS, sizeof(string::size_type));
	buf = new char[nCS + 1];
	f.read(buf, nCS + 1);
	msa->CS.assign(buf, nCS); /* read the null terminal */
	delete[] buf;
	/* read concatMSA */
	buf = new char[msa->numSeq * msa->csLen + 1];
	f.read(buf, msa->numSeq * msa->csLen + 1);
	msa->concatMSA.assign(buf, msa->numSeq * msa->csLen);
	delete[] buf;
	/* Read counts */
	msa->resCount = new unsigned*[msa->csLen];
	for(unsigned j = 0; j < msa->csLen; ++j) {
		msa->resCount[j] = new unsigned[msa->abc->getSize()];
		f.read((char*) msa->resCount[j], sizeof(unsigned) * msa->abc->getSize());
	}
	msa->gapCount = new unsigned[msa->csLen];
	f.read((char*) msa->gapCount, sizeof(unsigned) * msa->csLen);
	return msa;
}

} /* namespace EGriceLab */
