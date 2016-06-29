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
//#include "Stats.h"
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

    int32_t start = 1;
    int32_t end = concatLen;
	/* while there are possible occs and pattern not done */
    for (string::const_reverse_iterator rit = pattern.rbegin(); rit != pattern.rend() && start <= end; ++rit) {
      int8_t c = abc->encode(*rit) + 1; /* map pattern to our alphabet */
      start = C[c] + bwt->rank(c, start - 1); /* LF Mapping */
      end = C[c] + bwt->rank(c, end) - 1; /* LF Mapping */
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
      start = C[c] + bwt->rank(c, start - 1); /* LF Mapping */
      end = C[c] + bwt->rank(c, end) - 1; /* LF Mapping */
    }

    for(int32_t i = start; i <= end; ++i) {
    	uint32_t concatStart = accessSA(i);
    	int32_t csStart = concat2CS[concatStart];
    	int32_t csEnd = concat2CS[concatStart + pattern.length() - 1];
    	locs.push_back(CSLoc(csStart, csEnd, extractCS(concatStart, pattern)));
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
		start = C[c] + bwt->rank(c, start - 1); /* LF Mapping */
		end = C[c] + bwt->rank(c, end) - 1; /* LF Mapping */
		//cerr << "start:" << start << " end:" << end << endl;
	}

	if(start <= end) {
		uint32_t concatStart = accessSA(start);
    	int32_t csStart = concat2CS[concatStart];
    	int32_t csEnd = concat2CS[concatStart + pattern.length() - 1];
		return CSLoc(csStart, csEnd, extractCS(concatStart, pattern));
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
    	start = C[c] + bwt->rank(c, start - 1); /* LF Mapping */
    	end = C[c] + bwt->rank(c, end) - 1; /* LF Mapping */
    }
    if(start <= end) {
    	int32_t i = start + rand() % (end - start + 1);
    	uint32_t concatStart = accessSA(i); // sample a random position
    	int32_t csStart = concat2CS[concatStart];
    	int32_t csEnd = concat2CS[concatStart + pattern.length() - 1];
    	return CSLoc(csStart, csEnd, extractCS(concatStart, pattern));
    }
    else
    	return CSLoc();
}

bool CSFMIndex::save(ofstream& f) const {
	/* write alphabet name */
	unsigned nAlphabet = abc->getName().length();
	f.write((char*) &nAlphabet, sizeof(unsigned));
	f.write(abc->getName().c_str(), nAlphabet + 1); /* write the null terminal */

	/* write gap char */
	f.write(&gapCh, sizeof(char));

	/* write sizes */
	f.write((char*) &csLen, sizeof(uint16_t));
	f.write((char*) &concatLen, sizeof(int32_t));

	/* write arrays and objects */
	f.write((char*) C, (UINT8_MAX + 1) * sizeof(int32_t));

	f.write((char*) csSeq.c_str(), csLen + 2); /* write the null terminal */
	f.write((char*) csIdentity, (csLen + 1) * sizeof(double));
	f.write((char*) concat2CS, (concatLen + 1) * sizeof(uint16_t));
	f.write((char*) saSampled, (concatLen / SA_SAMPLE_RATE) * sizeof(uint32_t));

	saIdx->save(f);
	bwt->save(f);

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

	/* read gap char */
	f.read(&idx->gapCh, sizeof(char));

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

    idx->concat2CS = new uint16_t[idx->concatLen + 1];
	f.read((char*) idx->concat2CS, (idx->concatLen + 1) * sizeof(uint16_t));

	idx->saSampled = new uint32_t[idx->concatLen / SA_SAMPLE_RATE + 1];
	f.read((char*) idx->saSampled, (idx->concatLen / SA_SAMPLE_RATE) * sizeof(uint32_t));

	idx->saIdx = BitSequenceRRR::load(f); /* use RRR implementation */
	idx->bwt = WaveletTreeNoptrs::load(f);
	return idx;
}

CSFMIndex::CSFMIndex() : abc(NULL), gapCh('\0'), csLen(0),
		concatLen(0), C() /* zero-initiation */, csIdentity(), /* zero-initiation */
		concat2CS(NULL), saSampled(NULL), saIdx(NULL), bwt(NULL) {
}

CSFMIndex::~CSFMIndex() {
	//delete[] concatSeq;
	delete[] csIdentity;
	delete[] concat2CS;
	delete[] saSampled;
	delete saIdx;
	delete bwt;
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
	csFM->concatLen = msa->getMSANonGapLen() + msa->getNumSeq(); /* including one seperator per seq */
	csFM->csSeq = ' ' + msa->getCS(); /* dummy position 0 w/ white-space */
	csFM->csIdentity = new double[csFM->csLen + 1];
	for(unsigned j = 0; j < csFM->csLen; ++j)
		csFM->csIdentity[j + 1] = msa->identityAt(j);

	const int32_t N = csFM->concatLen + 1;

	/* construct the concatSeq and aux data */
	uint8_t* concatSeq = new uint8_t[N]; /* null terminated encoded string */
	csFM->concat2CS = new uint16_t[N](); /* zero-initiate 1-based concatSeq pos to CS pos, 0 for gap pos on CS */

	string::size_type shift = 0;
	for(unsigned i = 0; i < msa->getNumSeq(); ++i) {
		for(unsigned j = 0; j < csFM->csLen; ++j) {
			char c = msa->residualAt(i, j);
			if(!csFM->abc->isGap(c)) {
				int8_t k = csFM->abc->encode(::toupper(c)) + 1; /* encode to 1..alphabet-size range */
				csFM->C[k]++; // count alphabet frequency
				concatSeq[shift] = k; /* always store upper-case characters */
				csFM->concat2CS[shift] = j + 1; /* 1-based consensus position */
				shift++;
			}
		}
		csFM->C[sepCh]++; // count the separator
		concatSeq[shift] = sepCh; // add a separator
		csFM->concat2CS[shift] = 0; // separator point to gap
		shift++;
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

    /* construct SA */
    saidx_t errn;
    int32_t* SA = new int32_t[N];
	errn = divsufsort(concatSeq, SA, N);
	if(errn != 0)
		throw runtime_error("Error: Cannot build suffix-array on forward concatenated seq");

    /* construct the saSampled and saIdx */
	csFM->saSampled = new uint32_t[N / SA_SAMPLE_RATE + 1]();
	uint32_t* saHead = csFM->saSampled;
	BitString B(N); /* a temp BitString for building saIdx */
	for(uint32_t i = 0; i < N; ++i)
		if(SA[i] % SA_SAMPLE_RATE == 0) {
			*saHead++ = SA[i];
			B.setBit(i);
		}
	//cerr << "shift:" << saHead - csFM->saSampled << endl;

    csFM->saIdx = new BitSequenceRRR(B, RRR_SAMPLE_RATE);

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

    csFM->bwt = new WaveletTreeNoptrs((uint32_t *) X_bwt, N,
    		sizeof(uint8_t) * 8, bsb, map, true); // free the X_bwt after use

	return csFM;
}

uint32_t CSFMIndex::accessSA(uint32_t i) const {
	int32_t dist = 0;
	while(!saIdx->access(i)) {
		uint8_t c = bwt->access(i);
		i = C[c] + bwt->rank(c, i) - 1; // backward LF-mapping
		dist++;
	}
	return saSampled[saIdx->rank1(i) - 1] + dist;
}

string CSFMIndex::extractCS(int32_t start, const string& pattern) const {
	string csSeq;
	if(pattern.empty())
		return csSeq; // return empty CS
	assert(concat2CS[start] != 0 && concat2CS[start + pattern.length() - 1] != 0);
	for(int32_t i = start; i < start + pattern.length(); ++i) {
		if(i > start && concat2CS[i] - concat2CS[i - 1] > 1) // there are gaps between this two location
			csSeq.append(concat2CS[i] - concat2CS[i - 1] - 1, gapCh);
		csSeq.push_back(pattern[i - start]);
	}
	return csSeq;
}

} /* namespace EGriceLab */

