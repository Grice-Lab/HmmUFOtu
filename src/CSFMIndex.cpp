/*
 * CSFMIndex.cpp
 *
 *  Created on: Nov 5, 2015
 *      Author: zhengqi
 */
#include <algorithm>
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
      start = C[c] + fwt_bwt->rank(c, start - 1); /* LF Mapping */
      end = C[c] + fwt_bwt->rank(c, end) - 1; /* LF Mapping */
    }
    return start <= end ? end - start + 1 : 0;
}

vector<CSLoc> CSFMIndex::locate(const string& pattern) const {
	vector<CSLoc> locs;
	if(pattern.empty())
		return locs; /* empty pattern matches to nothing */
    int32_t start = 1; /* starting range in M from p[m-1] */
    int32_t end = concatLen;
	/* while there are possible occs and pattern not done */
    for (string::const_reverse_iterator rit = pattern.rbegin(); rit != pattern.rend() && start <= end; ++rit) {
      int8_t c = abc->encode(*rit) + 1; /* map pattern to 1 .. size */
      start = C[c] + fwt_bwt->rank(c, start - 1); /* LF Mapping */
      end = C[c] + fwt_bwt->rank(c, end) - 1; /* LF Mapping */
    }

    for(int32_t i = start; i <= end; ++i) {
    	int32_t csStart = csSA[i - 1];
    	int32_t csEnd = csSA[i - 1 + pattern.length() - 1];
    	locs.push_back(CSLoc(csStart, csEnd, extractCS(i, pattern.length())));
    }
    return locs;
}

CSLoc CSFMIndex::locateFirst(const string& pattern) const {
	if(pattern.empty())
		return CSLoc(); /* empty pattern matches to nothing */
	int32_t start = 1;
	int32_t end = concatLen;
	/* while there are possible occs and pattern not done */
	for (string::const_reverse_iterator rit = pattern.rbegin(); rit != pattern.rend() && start <= end; ++rit) {
		int8_t c = abc->encode(*rit) + 1; /* map pattern to 1 .. size */
		start = C[c] + fwt_bwt->rank(c, start - 1); /* LF Mapping */
		end = C[c] + fwt_bwt->rank(c, end) - 1; /* LF Mapping */
		//cerr << "start:" << start << " end:" << end << endl;
	}

	if(start <= end) {
    	int32_t csStart = csSA[start - 1];
    	int32_t csEnd = csSA[start - 1 + pattern.length() - 1];
		return CSLoc(csStart, csEnd, extractCS(start, pattern.length()));
	}
	else
		return CSLoc();
}

CSLoc CSFMIndex::locateOne(const string& pattern) const {
	if(pattern.empty())
		return CSLoc(); /* empty pattern matches to nothing */
    int32_t start = 1;
    int32_t end = concatLen;
	/* while there are possible occs and pattern not done */
    for (string::const_reverse_iterator rit = pattern.rbegin(); rit != pattern.rend() && start <= end; ++rit) {
    	int8_t c = abc->encode(*rit) + 1; /* map pattern to 1 .. size */
    	start = C[c] + fwt_bwt->rank(c, start - 1); /* LF Mapping */
    	end = C[c] + fwt_bwt->rank(c, end) - 1; /* LF Mapping */
    	//cerr << "start:" << start << " end:" << end << endl;
    }
    if(start <= end) {
    	int32_t i = start + rand() % (end - start + 1); // sample a random position between start and end
    	//int32_t i = start;
    	//uint8_t c = RRR_bwt->access(i);
    	//int32_t j = C[c] + RRR_bwt->rank(c, i) - 1; /* LF-mapping */
    	//cerr << "i:" << i << " c:" << abc->decode(c) << " j:" << j << " csSA:" << csSA[j] << endl;
    	int32_t csStart = csSA[i - 1];
    	int32_t csEnd = csSA[i - 1 + pattern.length() - 1];
    	return CSLoc(csStart, csEnd, extractCS(i, pattern.length()));
    }
    else
    	return CSLoc();
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
	f.write((char*) csSA, (concatLen + 1) * sizeof(uint16_t));

	fwt_bwt->save(f);
	rev_bwt->save(f);

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

    idx->csSA = new uint16_t[idx->concatLen + 1];
	f.read((char*) idx->csSA, (idx->concatLen + 1) * sizeof(uint16_t));

	idx->fwt_bwt = WaveletTreeNoptrs::load(f);
	idx->rev_bwt = WaveletTreeNoptrs::load(f);

	return idx;
}

CSFMIndex::CSFMIndex() : abc(NULL), gapCh('\0'), csLen(0),
		concatLen(0), C() /* zero-initiation */, csIdentity(), /* zero-initiation */
		csSA(NULL), fwt_bwt(NULL), rev_bwt(NULL) {
}

CSFMIndex::~CSFMIndex() {
	//delete[] concatSeq;
	delete[] csIdentity;
	delete[] csSA;
	delete[] fwt_bwt;
	delete[] rev_bwt;
}

CSFMIndex* CSFMIndex::build(const MSA* msa) {
	if(!(msa->getCSLen() <= UINT16_MAX))
		throw invalid_argument("CSFMIndex cannot handle MSA with consensus length longer than " + UINT16_MAX);
	CSFMIndex* csFM = new CSFMIndex; // allocate an empty object
	//cerr << "FM initiated" << endl;
	/* set basic info */
	csFM->abc = msa->getAbc();
	csFM->gapCh = csFM->abc->getGap()[0]; // use the default gap character
	csFM->csLen = msa->getCSLen();
	csFM->concatLen = msa->getMSANonGapLen();
	csFM->csSeq = ' ' + msa->getCS(); /* dummy position 0 w/ white-space */
	csFM->csIdentity = new double[csFM->csLen + 1];
	for(unsigned j = 0; j < csFM->csLen; ++j)
		csFM->csIdentity[j + 1] = msa->identityAt(j);

	const int32_t N = csFM->concatLen + 1;

	/* construct the concatSeq and aux data */
	uint8_t* concatSeq = new uint8_t[N]; /* null terminated encoded string */
	uint16_t* concat2CS = new uint16_t[N](); /* zero-initiate 1-based concatSeq pos to CS pos, 0 for gap pos on CS */

	string::size_type shift = 0;
	for(unsigned i = 0; i < msa->getNumSeq(); ++i) {
		for(unsigned j = 0; j < csFM->csLen; ++j) {
			char c = msa->residualAt(i, j);
			if(!csFM->abc->isGap(c)) {
				int8_t k = csFM->abc->encode(::toupper(c)) + 1; /* encode to 1..alphabet-size range */
				csFM->C[k]++; // count alphabet frequency
				concatSeq[shift] = k; /* always store upper-case characters */
				concat2CS[shift] = j + 1; /* 1-based consensus position */
				shift++;
			}
		}
	}
	assert(shift == N - 1);
	concatSeq[shift] = '\0'; // null terminal
	csFM->C['\0']++; // count the null terminal

	/* construct cumulative counts */
    int32_t prev = csFM->C[0];
    int32_t tmp;
    csFM->C[0] = 0;
    for (int i = 1; i <= csFM->abc->getSize() + 1; ++i) {
      tmp = csFM->C[i];
      csFM->C[i] = csFM->C[i-1] + prev;
      prev = tmp;
    }

	/* construct the fwd index */
    /* construct SA */
    saidx_t errno;
    int32_t* SA = new int32_t[N];
	errno = divsufsort(concatSeq, SA, N);
	if(errno != 0)
		throw runtime_error("Error: Cannot build suffix-array on forward concatenated seq");

    /* construct primary index that map loc on SA/FM-indx to CS index */
    csFM->csSA = new uint16_t[N];
    for(int32_t i = 0; i < N; ++i)/* check each SA pos */
    	csFM->csSA[i] = concat2CS[SA[i]];

    delete[] concat2CS; // delete concat2CS once csSA is built

    /* construct BWT and index */
	uint8_t* X_bwt = new uint8_t[N];
	if(X_bwt == NULL)
		throw runtime_error("Error: Cannot allocate BWT string for concatSeq");
    for(int32_t i = 0; i < N; ++i)
        if(SA[i] == 0) // matches to the null
            X_bwt[i] = '\0'; // null terminal
        else X_bwt[i] = concatSeq[SA[i] - 1];

	/* construct RRR_compressed BWT */
    Mapper* map = new MapperNone(); /* smart ptr no delete necessary */
	BitSequenceBuilder* bsb = new BitSequenceBuilderRRR(RRR_SAMPLE_RATE); /* smart ptr no delete necessary */

    csFM->fwt_bwt = new WaveletTreeNoptrs((uint32_t *) X_bwt, N,
    		sizeof(uint8_t) * 8, bsb, map, false); // don't free the X_bwt yet, re-use later

	cerr << "Forward index built" << endl;

	/* construct the rev index */
	std::reverse(concatSeq, concatSeq + N - 1); // re-use the concatSeq memory

    /* construct SA */
	errno = divsufsort(concatSeq, SA, N); // re-use the SA
	if(errno != 0)
		throw runtime_error("Error: Cannot build suffix-array on reverse concatenated seq");

    /* construct BWT and index, re-use X-bwt */
    for(int32_t i = 0; i < N; ++i)
        if(SA[i] == 0) // matches to the null
            X_bwt[i] = '\0'; // null terminal
        else X_bwt[i] = concatSeq[SA[i] - 1];

	/* construct RRR_compressed BWT, re-use memory */
    csFM->rev_bwt = new WaveletTreeNoptrs((uint32_t *) X_bwt, N,
    		sizeof(uint8_t) * 8, bsb, map, true);

	cerr << "Reverse index built" << endl;
    //delete[] X_bwt; /* X_bwt already deleted during WaveletTreeNoptrs construction */

	delete[] concatSeq;
    delete[] SA;

	return csFM;
}

string CSFMIndex::extractCS(int32_t start, int32_t len) const {
	string csSeq;
	if(len == 0)
		return csSeq; // return empty CS
	for(int32_t i = start; i < start + len; ++i) {
		if(i > 0 && csSA[i] - csSA[i - 1] > 1) // there are gaps between this two location
			csSeq.append(csSA[i] - csSA[i - 1] - 1, gapCh);
		csSeq.push_back(abc->decode(fwt_bwt->access(i)));
	}
	return csSeq;
}

} /* namespace EGriceLab */

