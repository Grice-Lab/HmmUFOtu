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

set<unsigned> CSFMIndex::locateIndex(const string& pattern) const {
    set<unsigned> idx;
	if(pattern.empty())
		return idx; /* empty pattern matches to nothing */
    int32_t start = 1;
    int32_t end = concatLen;
	/* while there are possible occs and pattern not done */
    for (string::const_reverse_iterator rit = pattern.rbegin(); rit != pattern.rend() && start <= end; ++rit) {
    	int8_t c = abc->encode(*rit) + 1; /* map pattern to 1 .. size */
    	start = C[c] + bwt->rank(c, start - 1); /* LF Mapping */
    	end = C[c] + bwt->rank(c, end) - 1; /* LF Mapping */
    }

    for(int32_t i = start; i <= end; ++i) {
    	uint32_t k = accessSA(i); /* concatStart */
    	idx.insert(k / (csLen + 1));
    }

    return idx;
}

ofstream& CSFMIndex::save(ofstream& out) const {
	/* write alphabet name */
	unsigned nAlphabet = abc->getName().length();
	out.write((char*) &nAlphabet, sizeof(unsigned));
	out.write(abc->getName().c_str(), nAlphabet + 1); /* write the null terminal */

	/* write gap char */
	out.write(&gapCh, sizeof(char));

	/* write sizes */
	out.write((char*) &csLen, sizeof(uint16_t));
	out.write((char*) &concatLen, sizeof(int32_t));

	/* write arrays and objects */
	out.write((char*) C, (UINT8_MAX + 1) * sizeof(int32_t));

	out.write((char*) csSeq.c_str(), csLen + 2); /* write the null terminal */
	out.write((char*) csIdentity, (csLen + 1) * sizeof(double));
	out.write((char*) concat2CS, (concatLen + 1) * sizeof(uint16_t));
	out.write((char*) saSampled, (concatLen / SA_SAMPLE_RATE) * sizeof(uint32_t));

	saIdx->save(out);
	bwt->save(out);

	return out;
}

ifstream& CSFMIndex::load(ifstream& in) {
	clear(); /* clear old data, if any */
	/* read alphabet by name */
	char* buf = NULL; /* buf for reading */
	unsigned nAlphabet;
	in.read((char *) &nAlphabet, sizeof(unsigned));
	buf = new char[nAlphabet + 1]; /* include the null terminal */
	in.read(buf, nAlphabet + 1);
	abc = AlphabetFactory::getAlphabetByName(buf);
	delete[] buf;

	/* read gap char */
	in.read(&gapCh, sizeof(char));

	/* read sizes */
	in.read((char*) &csLen, sizeof(uint16_t));
	in.read((char*) &concatLen, sizeof(int32_t));

	/* read arrays and objects */
	in.read((char*) C, (UINT8_MAX + 1) * sizeof(int32_t));

	buf = new char[csLen + 2];
	in.read((char*) buf, csLen + 2); /* read the null terminal */
	csSeq.assign(buf, csLen + 1); /* use assign to prevent memory leak */
	delete[] buf;

	csIdentity = new double[csLen + 1];
	in.read((char*) csIdentity, (csLen + 1) * sizeof(double));

    concat2CS = new uint16_t[concatLen + 1];
	in.read((char*) concat2CS, (concatLen + 1) * sizeof(uint16_t));

	saSampled = new uint32_t[concatLen / SA_SAMPLE_RATE + 1];
	in.read((char*) saSampled, (concatLen / SA_SAMPLE_RATE) * sizeof(uint32_t));

	saIdx = BitSequenceRRR::load(in); /* use RRR implementation */
	bwt = WaveletTreeNoptrs::load(in);
	return in;
}

CSFMIndex& CSFMIndex::build(const MSA& msa) {
	if(!(msa.getCSLen() <= UINT16_MAX)) {
		throw runtime_error("CSFMIndex cannot handle MSA with consensus length longer than " + UINT16_MAX);
		return *this;
	}
	clear(); /* clear old data, if any */

	/* construct basic information */
	buildBasic(msa);
	/* construct concatSeq related info */
	uint8_t* concatSeq = buildConcatSeq(msa);

    /* construct SA and BWT */
    buildBWT(concatSeq);

    /* free temporary memories */
    delete[] concatSeq;
	return *this;
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

void CSFMIndex::buildBasic(const MSA& msa) {
	abc = msa.getAbc();
	gapCh = abc->getGap().front(); // use the default gap character
	csLen = msa.getCSLen();
	concatLen = msa.getMSANonGapLen() + msa.getNumSeq(); /* including one seperator per seq */
	csSeq = ' ' + msa.getCS(); /* dummy position 0 w/ white-space */
	csIdentity = new double[csLen + 1];
	csIdentity[0] = 0; /* dummy value */
	for(unsigned j = 0; j < csLen; ++j)
		csIdentity[j + 1] = msa.identityAt(j);
}

uint8_t* CSFMIndex::buildConcatSeq(const MSA& msa) {
	const int32_t N = concatLen + 1;
	/* construct the concatSeq update concat2CS index */
	uint8_t* concatSeq = new uint8_t[N]; /* null terminated encoded string */
	concat2CS = new uint16_t[N](); /* zero-initiate 1-based concatSeq pos to CS pos, 0 for gap pos on CS */

	string::size_type shift = 0;
	for(unsigned i = 0; i < msa.getNumSeq(); ++i) {
		for(unsigned j = 0; j < csLen; ++j) {
			char c = msa.residualAt(i, j);
			if(!abc->isGap(c)) {
				int8_t k = abc->encode(::toupper(c)) + 1; /* encode to 1..alphabet-size range */
				C[k]++; // count alphabet frequency
				concatSeq[shift] = k; /* always store upper-case characters */
				concat2CS[shift] = j + 1; /* 1-based consensus position */
				shift++;
			}
		}
		C[sepCh]++; // count the separator
		concatSeq[shift] = sepCh; // add a separator at the end of each seq
		concat2CS[shift] = 0; // separator point to gap
		shift++;
	}
	assert(shift == N - 1);
	concatSeq[shift] = '\0'; // add null terminal
	C['\0']++; // count the null terminal

	/* construct cumulative counts */
    int32_t prev = C[0];
    int32_t tmp;
    C[0] = 0;
    for (int i = 1; i <= abc->getSize() + 1; ++i) {
      tmp = C[i];
      C[i] = C[i-1] + prev;
      prev = tmp;
    }

	return concatSeq;
}

void CSFMIndex::buildBWT(const uint8_t* concatSeq) {
    /* construct SA */
    saidx_t errn;
    const int32_t N = concatLen + 1;
    int32_t* SA = new int32_t[N];

	errn = divsufsort(concatSeq, SA, N);
	if(errn != 0)
		throw runtime_error("Error: Cannot build suffix-array on forward concatenated seq");

    /* construct the saSampled and saIdx */
	saSampled = new uint32_t[N / SA_SAMPLE_RATE + 1]() /* zero-initiation */;
	uint32_t* saHead = saSampled;
	BitString B(N); /* a temp BitString for building saIdx */
	for(uint32_t i = 0; i < N; ++i)
		if(SA[i] % SA_SAMPLE_RATE == 0) {
			*saHead++ = SA[i];
			B.setBit(i);
		}
	//cerr << "shift:" << saHead - csFM->saSampled << endl;

    saIdx = new BitSequenceRRR(B, RRR_SAMPLE_RATE); /* use RRR implementation */

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
	BitSequenceBuilder* bsb = new BitSequenceBuilderRRR(RRR_SAMPLE_RATE); /* bsb is a smart ptr no delete necessary */

    bwt = new WaveletTreeNoptrs((uint32_t *) X_bwt, N,
    		sizeof(uint8_t) * 8, bsb, map, true); // free the X_bwt after use

    /* free temporary memories */
    delete[] SA;
}


} /* namespace EGriceLab */

