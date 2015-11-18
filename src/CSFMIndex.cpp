/*
 * CSFMIndex.cpp
 *
 *  Created on: Nov 5, 2015
 *      Author: zhengqi
 */
#include <stdexcept>
#include <stdint.h>
#include <cctype>
#include <cstdlib>
#include "CSFMIndex.h"
#include "Stats.h"
#include "BitSequenceBuilder.h"
#include "BitSequenceBuilderRRR.h"
#include "Mapper.h"
#include "MapperNone.h"

using namespace std;
using namespace cds_static;
namespace EGriceLab {

int32_t CSFMIndex::count(const string& pattern) const {
	int32_t m = pattern.length();
	if(m == 0)
		return 0; /* empty pattern matches to nothing */
    uint8_t c = abc->encode(pattern[m-1]); /* map pattern to our alphabet */
    int32_t i = m - 1;

    int32_t start = C[c]; /* starting range in M from p[m-1] */
    int32_t end = C[c+1]-1;
	/* while there are possible occs and pattern not done */
    while (start <= end && i >= 1) {
      c = abc->encode(pattern[--i]); /* map pattern to our alphabet */
      start = C[c] + RRR_bwt->rank(c, start - 1); /* LF Mapping */
      end = C[c] + RRR_bwt->rank(c, end) - 1; /* LF Mapping */
    }
    return start <= end ? end - start + 1 : 0;
}

vector<uint16_t> CSFMIndex::locate(const string& pattern) const {
	vector<uint16_t> locs;
	int32_t m = pattern.length();
	if(m == 0)
		return locs; /* empty pattern matches to nothing */
    uint8_t c = abc->encode(pattern[m-1]); /* map pattern to our alphabet */
    int32_t i = m - 1;

    int32_t start = C[c]; /* starting range in M from p[m-1] */
    int32_t end = C[c+1]-1;
	/* while there are possible occs and pattern not done */
    while (start <= end && i >= 1) {
      c = abc->encode(pattern[--i]); /* map pattern to our alphabet */
      start = C[c] + RRR_bwt->rank(c, start - 1); /* LF Mapping */
      end = C[c] + RRR_bwt->rank(c, end) - 1; /* LF Mapping */
    }

    for(int32_t i = start; i <= end; ++i) {
    	c = RRR_bwt->access(i);
    	int32_t j = C[c] + RRR_bwt->rank(c, i) - 1; /* LF-mapping */
    	locs.push_back(csSA[j]);
    }
    return locs;
}

uint16_t CSFMIndex::locateFirst(const string& pattern) const {
	int32_t m = pattern.length();
	if(m == 0)
		return 0; /* empty pattern matches to nothing */
    uint8_t c = abc->encode(pattern[m-1]); /* map pattern to our alphabet */
    int32_t i = m - 1;

    int32_t start = C[c]; /* starting range in M from p[m-1] */
    int32_t end = C[c+1]-1;
	/* while there are possible occs and pattern not done */
    while (start <= end && i >= 1) {
      c = abc->encode(pattern[--i]); /* map pattern to our alphabet */
      start = C[c] + RRR_bwt->rank(c, start - 1); /* LF Mapping */
      end = C[c] + RRR_bwt->rank(c, end) - 1; /* LF Mapping */
    }

    if(start <= end) {
    	c = RRR_bwt->access(start);
    	int32_t j = C[c] + RRR_bwt->rank(c, start) - 1; /* LF-mapping */
    	return csSA[j];
    }
    else
    	return 0;
}

uint16_t CSFMIndex::locateOne(const string& pattern) const {
	int32_t m = pattern.length();
	if(m == 0)
		return concatLen; /* empty pattern matches to everywhere */
    uint8_t c = abc->encode(pattern[m-1]); /* map pattern to our alphabet */
    int32_t i = m - 1;

    int32_t start = C[c]; /* starting range in M from p[m-1] */
    int32_t end = C[c+1]-1;
	/* while there are possible occs and pattern not done */
    while (start <= end && i >= 1) {
      c = abc->encode(pattern[--i]); /* map pattern to our alphabet */
      start = C[c] + RRR_bwt->rank(c, start - 1); /* LF Mapping */
      end = C[c] + RRR_bwt->rank(c, end) - 1; /* LF Mapping */
    }

    if(start <= end) {
    	i = start + rand() % (end - start);
    	c = RRR_bwt->access(i);
    	int32_t j = C[c] + RRR_bwt->rank(c, i) - 1; /* LF-mapping */
    	return csSA[j];
    }
    else
    	return 0;
}

bool CSFMIndex::save(ofstream& f) const {
	/* write alphabet name */
	unsigned nAlphabet = abc->getName().length();
	f.write((char*) &nAlphabet, sizeof(unsigned));
	f.write(abc->getName().c_str(), nAlphabet + 1); /* write the null terminal */

	/* write sizes */
	f.write((char*) &csLen, sizeof(uint16_t));
	f.write((char*) &concatLen, sizeof(int32_t));

	/* write arrays and objects */
	f.write((char*) C, (UINT8_MAX + 1) * sizeof(int32_t));

	f.write((char*) csSeq.c_str(), csLen + 2); /* write the null terminal */
	f.write((char*) csIdentity, (csLen + 1) * sizeof(double));
	f.write((char*) csSA, concatLen * sizeof(uint16_t));

	RRR_bwt->save(f);

	return f.good();
}

CSFMIndex* CSFMIndex::load(ifstream& f) {
	CSFMIndex* idx = new CSFMIndex; /* Allocate a new object */
	/* read alphabet by name */
	char* buf = NULL; /* buf for reading */
	unsigned nAlphabet;
	f.read((char *) &nAlphabet, sizeof(unsigned));
	buf = new char[nAlphabet + 1]; /* include the null terminal */
	f.read(buf, nAlphabet + 1);
	idx->abc = SeqCommons::getAlphabetByName(buf);
	delete[] buf;

	/* read sizes */
	f.read((char*) &idx->csLen, sizeof(uint16_t));
	f.read((char*) &idx->concatLen, sizeof(int32_t));

	/* read arrays and objects */
	f.read((char*) idx->C, (UINT8_MAX + 1) * sizeof(int32_t));

	buf = new char[idx->csLen + 2];
	f.read((char*) buf, idx->csLen + 2); /* read the null terminal */
	idx->csSeq.assign(buf, idx->csLen + 2); /* use assign to prevent memory leak */
	delete[] buf;

	idx->csIdentity = new double[idx->csLen + 1];
	f.read((char*) idx->csIdentity, (idx->csLen + 1) * sizeof(double));

    idx->csSA = new uint16_t[idx->concatLen];
	f.read((char*) idx->csSA, idx->concatLen * sizeof(uint16_t));

	idx->RRR_bwt = WaveletTreeNoptrs::load(f);

	return idx;
}

CSFMIndex::CSFMIndex() : abc(NULL), csLen(0),
		concatLen(0), C() /* zero-initiation */, csIdentity(), /* zero-initiation */
		csSA(NULL), RRR_bwt(NULL) {
}

CSFMIndex::~CSFMIndex() {
	//delete[] concatSeq;
	delete[] csIdentity;
	delete[] csSA;
	delete[] RRR_bwt;
}

CSFMIndex* CSFMIndex::build(const MSA* msa) {
	if(!(msa->getCSLen() <= UINT16_MAX + 1))
		throw invalid_argument("CSFMIndex cannot handle MSA with consensus length longer than " + UINT16_MAX + 1);
	CSFMIndex* csFM = new CSFMIndex; // allocate an empty object
	//cerr << "FM initiated" << endl;
	/* set basic info */
	csFM->abc = msa->getAbc();
	csFM->csLen = msa->getCSLen();
	csFM->concatLen = msa->getMSANonGapLen();
	csFM->csSeq = ' ' + msa->getCS(); /* dummy position 0 */
	csFM->csIdentity = new double[csFM->csLen + 1];
	for(unsigned j = 0; j < csFM->csLen; ++j)
		csFM->csIdentity[j + 1] = msa->identityAt(j);

	//printVerboseInfo("basic info built");

	/* construct the concatSeq and aux data */
	uint8_t* concatSeq = new uint8_t[csFM->concatLen]; /* non-null terminated string */
	uint16_t* concat2CS = new uint16_t[csFM->concatLen](); /* zero-initiate 1-based concatSeq pos to CS pos, 0 for gap pos on CS */

	string::size_type shift = 0;
	for(unsigned i = 0; i < msa->getNumSeq(); ++i) {
		for(unsigned j = 0; j < csFM->csLen; ++j) {
			char c = msa->residualAt(i, j);
			if(!csFM->abc->isGap(c)) {
				concatSeq[shift] = csFM->abc->encode(::toupper(c)); /* always store upper-case characters */
				concat2CS[shift] = j + 1; /* 1-based consensus position */
				shift++;
			}
		}
	}
	assert(shift == csFM->concatLen);

	/* construct BWT and index */
	//uint8_t* X_bwt = new uint8_t[concatLen];
	int32_t* bwtIdx = new int32_t[csFM->concatLen]; /* index mapping pos on Idx to pos on cocnatSeq */
	uint8_t* X_bwt = concatSeq; /* make a new alias */
	saidx_t errno = bw_transform(concatSeq, X_bwt, NULL, csFM->concatLen, bwtIdx); /* in-place BWT transform */
	if(errno != 0)
		throw runtime_error("Error: Cannot build BWT and index on concatenated seq");

	//cerr << "BWT transformation done" << endl;

	/* construct RRR_compressed BWT */
    Mapper* map = new MapperNone(); /* smart ptr no delete necessary */
	BitSequenceBuilder* bsb = new BitSequenceBuilderRRR(RRR_SAMPLE_RATE); /* smart ptr no delete necessary */

    csFM->RRR_bwt = new WaveletTreeNoptrs((uint *) X_bwt, csFM->concatLen / (sizeof(uint) / sizeof(uint8_t)),
    		sizeof(uint8_t) * sizeof(char), bsb, map, true);

	//cerr << "RRR transformation done" << endl;

    /* Transfer the primary index to CS index */
    csFM->csSA = new uint16_t[csFM->concatLen];
    for(int32_t i = 0; i < csFM->concatLen; ++i) /* check each SA pos */
    	csFM->csSA[i] = concat2CS[bwtIdx[i]];
    //cerr << "csSA built" << endl;
    /* Delete local arrays */
    delete[] concat2CS;
    //cerr << "concat2CS deleted" << endl;
    //delete[] X_bwt; /* X_bwt already deleted during WaveletTreeNoptrs construction */
    delete[] bwtIdx;

	return csFM;
}

} /* namespace EGriceLab */

