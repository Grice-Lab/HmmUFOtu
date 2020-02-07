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
 * HmmUFOtu_main.cpp
 *  source file for HmmUFOtu core algorithms
 *  Created on: Jul 10, 2017
 *      Author: zhengqi
 */

#include <Eigen/Dense>
#include <cassert>
#include <algorithm>
#include "HmmUFOtu_main.h"
#include "StringUtils.h"

using namespace std;
using namespace Eigen;

namespace EGriceLab {
namespace HmmUFOtu {

BandedHMMP7::HmmAlignment alignSeq(const BandedHMMP7& hmm, const CSFMIndex& csfm, const PrimarySeq& read,
		int seedLen, int seedRegion, BandedHMMP7::align_mode mode) {
	const DegenAlphabet* abc = hmm.getNuclAbc();
	const int K = hmm.getProfileSize();
	const int L = hmm.getCSLen();
	const int N = read.length();

	BandedHMMP7::ViterbiScores seqVscore(K, N); // construct an empty reusable score
	vector<BandedHMMP7::ViterbiAlignPath> seqVpaths; // construct an empty list of VPaths
	BandedHMMP7::ViterbiAlignTrace seqVtrace; // construct an empty VTrace

	int regionLen = seedRegion < read.length() ? seedRegion : read.length(); /* search region */
	/* find seed in 5' */
	for(int seedFrom = 0; seedFrom + seedLen - 1 < regionLen; ++seedFrom) {
		int seedTo = seedFrom + seedLen - 1;
		PrimarySeq seed(abc, read.getId(), read.subseq(seedFrom, seedLen));
		const CSLoc& loc = csfm.locateOne(seed.getSeq());
		if(loc.isValid()) /* a read seed located */ {
//			cerr << "using 5' seed seedFrom: " << seedFrom << " seedTo: " << seedTo << endl;
//			cerr << "Using 5' seed: " << seed.getSeq() << endl;
//			fprintf(stderr, "start:%d end:%d from:%d to:%d  CSLen:%d CS:%s\n", loc.start, loc.end, seedFrom + 1, seedFrom + seedLen, loc.CS.length(), loc.CS.c_str());
			const BandedHMMP7::ViterbiAlignPath& vpath = hmm.buildAlignPath(loc, seedFrom + 1, seedTo + 1);
			if(vpath.isValid()) {
				seqVpaths.push_back(vpath); /* seed_from and seed_to are 1-based */
				break; /* only one 5'-seed necessary */
			}
		}
	}
	/* find seed in 3', if requested */
	if(mode == BandedHMMP7::GLOBAL && (seqVpaths.empty() || read.length() >= 2 * regionLen)) {
		for(int seedTo = read.length() - 1; seedTo - seedLen + 1 >= (int) read.length() - regionLen; --seedTo) {
			int seedFrom = seedTo - seedLen + 1;
			PrimarySeq seed(abc, read.getId(), read.subseq(seedFrom, seedLen));
			const CSLoc& loc = csfm.locateOne(seed.getSeq());
			if(loc.isValid()) { /* a read seed located */
//				cerr << "using 3' seed seedFrom: " << seedFrom << " seedTo: " << seedTo << endl;
//				cerr << "Using 3' seed: " << seed.getSeq() << endl;
//				fprintf(stderr, "start:%d end:%d from:%d to:%d  CSLen:%d CS:%s\n", loc.start, loc.end, seedTo - seedLen + 2, seedTo + 1, loc.CS.length(), loc.CS.c_str());
				const BandedHMMP7::ViterbiAlignPath& vpath = hmm.buildAlignPath(loc, seedFrom + 1, seedTo + 1);
				if(vpath.isValid()) {
					seqVpaths.push_back(vpath); /* seed_from and seed_to are 1-based */
					break; /* only one 3'-seed necessary */
				}
			}
		}
	}

	/* banded HMM align */
	if(!seqVpaths.empty()) { /* use banded Viterbi algorithm */
		hmm.calcViterbiScores(read, seqVscore, seqVpaths);
		if(seqVscore.S.minCoeff() == inf) { /* banded version failed */
			debugLog << "Banded HMM algorithm didn't find a potential Viterbi path, returning to regular HMM" << endl;
			seqVscore.reset();
			hmm.calcViterbiScores(read, seqVscore);
		}
	}
	else
		hmm.calcViterbiScores(read, seqVscore); /* use original Viterbi algorithm */

	/* build VTrace */
	hmm.buildViterbiTrace(seqVscore, seqVtrace);

	assert(seqVtrace.minScore != inf);

	/* get aligned seq */
	return hmm.buildGlobalAlign(read, seqVscore, seqVtrace);
}

BandedHMMP7::HmmAlignment alignSeq(const BandedHMMP7& hmm, const PrimarySeq& read) {
	const DegenAlphabet* abc = hmm.getNuclAbc();
	const int K = hmm.getProfileSize();
	const int L = hmm.getCSLen();
	const int N = read.length();

	BandedHMMP7::ViterbiScores seqVscore(K, N); // construct an empty reusable score
	BandedHMMP7::ViterbiAlignTrace seqVtrace; // construct an empty VTrace

	/* traditional HMM align */
	hmm.calcViterbiScores(read, seqVscore); /* use original Viterbi algorithm */

	/* build VTrace */
	hmm.buildViterbiTrace(seqVscore, seqVtrace);

	assert(seqVtrace.minScore != inf);
	/* get aligned seq */
	return hmm.buildGlobalAlign(read, seqVscore, seqVtrace);
}

vector<PTUnrooted::PTLoc> getSeed(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, double maxDiff, double maxHeight) {
	vector<PTUnrooted::PTLoc> locs; /* candidate locations */
	/* get potential placement locations based on pDist to observed or inferred sequences */
	for(vector<PTUnrooted::PTUNodePtr>::size_type i = 0; i < ptu.numNodes(); ++i) {
		PTUnrooted::PTUNodePtr node = ptu.getNode(i);
		if(!node->isRoot() && ptu.getHeight(node) <= maxHeight) {
			double pDist = SeqUtils::pDist(node->getSeq(), seq, start, end);
			locs.push_back(PTUnrooted::PTLoc(start, end, node->getId(), pDist));
		}
	}
	assert(!locs.empty());
	std::sort(locs.begin(), locs.end()); /* sort by p-Dist */
	/* remove bad seed, if necessary */
	double bestDist = locs[0].dist;
	double worstDist = locs[locs.size() - 1].dist;
	if(worstDist < bestDist + maxDiff) { /* need filtering */
		vector<PTUnrooted::PTLoc>::iterator goodSeed;
		for(goodSeed = locs.begin(); goodSeed != locs.end(); ++goodSeed) {
			if(goodSeed->dist - bestDist > maxDiff)
				break;
		}
		locs.erase(goodSeed, locs.end()); /* remove too bad placements */
	}
	return locs;
}

vector<PTUnrooted::PTPlacement> estimateSeq(const PTUnrooted& ptu, const DigitalSeq& seq,
		const vector<PTUnrooted::PTLoc>& locs, const string& method) {
	vector<PTUnrooted::PTPlacement> places;
	for(vector<PTUnrooted::PTLoc>::const_iterator loc = locs.begin(); loc != locs.end(); ++loc)
		places.push_back(ptu.estimateSeq(seq, *loc, method));
	return places;
}

vector<PTUnrooted::PTPlacement>& filterPlacements(vector<PTUnrooted::PTPlacement>& places, double maxError) {
	assert(!places.empty() && maxError >= 0);
	std::sort(places.rbegin(), places.rend(), compareByLoglik); /* sort places decently by estimated loglik */
	double bestEstLoglik = places[0].loglik;
	vector<PTUnrooted::PTPlacement>::iterator goodPlace;
	for(goodPlace = places.begin(); goodPlace != places.end(); ++goodPlace) {
		if(bestEstLoglik - goodPlace->loglik > maxError)
			break;
	}
	places.erase(goodPlace, places.end()); /* remove bad placements */
	return places;
}

vector<PTUnrooted::PTPlacement>& placeSeq(const PTUnrooted& ptu, const DigitalSeq& seq,
		vector<PTUnrooted::PTPlacement>& places, double maxHeight) {
	for(vector<PTUnrooted::PTPlacement>::iterator place = places.begin(); place != places.end(); ++place)
		ptu.placeSeq(seq, *place, maxHeight);
	return places;
}

void calcQValues(vector<PTUnrooted::PTPlacement>& places, PTUnrooted::PRIOR_TYPE type) {
	if(places.empty())
		return;

	/* explore all placements */
	VectorXd ppPlace(places.size()); /* posterior logP at placement */
	map<string, double> ppTaxon; /* posterior logP at taxon */
	double ppTaxNorm = infV; /* log(0) */

	VectorXd::Index i = 0;
	for(vector<PTUnrooted::PTPlacement>::const_iterator placement = places.begin(); placement != places.end(); ++placement) {
		double p = placement->loglik + placement->logPriorPr(type);
		ppPlace(i++) = p;
		string taxonomy = placement->getTaxonName();
		if(ppTaxon.find(taxonomy) == ppTaxon.end())
			ppTaxon[taxonomy] = p;
		else
			ppTaxon[taxonomy] = EGriceLab::Math::add_scaled(ppTaxon[taxonomy], p);
		ppTaxNorm = EGriceLab::Math::add_scaled(ppTaxNorm, p);
	}
	/* scale and normalize llPlace */
	VectorXd p = (ppPlace.array() - ppPlace.maxCoeff()).exp();
	p /= p.sum();
	/* calculate qPlace */
	for(vector<PTUnrooted::PTPlacement>::size_type i = 0; i < places.size(); ++i) {
		double q = EGriceLab::Math::p2q(1 - p(i));
		places[i].qPlace = q > PTUnrooted::PTPlacement::MAX_Q ? PTUnrooted::PTPlacement::MAX_Q : q;
	}

	/* calculate qTaxonomy */
	for(vector<PTUnrooted::PTPlacement>::iterator placement = places.begin(); placement != places.end(); ++placement) {
		double q = EGriceLab::Math::p2q(1 - ::exp(ppTaxon[placement->getTaxonName()] - ppTaxNorm));
		placement->qTaxon = q > PTUnrooted::PTPlacement::MAX_Q ? PTUnrooted::PTPlacement::MAX_Q : q;
	}
}

double alignIdentity(const DegenAlphabet* abc, const string& align, int start, int end) {
	assert(0 <= start && start <= end && end < align.size());
	int identity = 0;
	for(int i = start; i <= end; ++i)
		if(abc->isSymbol(align[i]))
			identity++;
	return static_cast<double> (identity) / (end - start + 1);
}

double hmmIdentity(const BandedHMMP7& hmm, const string& align, int start, int end) {
	assert(0 <= start && start <= end && end < align.size());
	int identity = 0;
	int nSite = 0;
	for(int i = start; i <= end; ++i) {
		if(hmm.getProfileLoc(i + 1) != 0) { /* a profile site */
			nSite++;
			if(hmm.getNuclAbc()->isSymbol(align[i]))
				identity++;
		}
	}
	return static_cast<double> (identity) / nSite;
}

JPlace::JPlace(int edgeID, string readName, double edgeLen, double ratio,
		double loglik, double annoDist, double q)
: edgeID(edgeID), readName(readName), likelihood(loglik), distal_length(edgeLen * ratio), proximal_length(edgeLen * (1.0 - ratio))
{
	pendant_length = ratio <= 0.5 ? annoDist - distal_length : annoDist - proximal_length;
	like_ratio = q >= MAX_Q ? 1 : Math::q2p(q);
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
