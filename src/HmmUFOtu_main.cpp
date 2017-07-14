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
#include "HmmUFOtu_main.h"
using namespace std;
using namespace Eigen;

namespace EGriceLab
{
const string PTPlacement::UNASSIGNED_TAXONNAME = "UNASSIGNED";
const double PTPlacement::UNASSIGNED_LOGLIK = EGriceLab::nan;
const string PTPlacement::UNASSIGNED_ID = "NULL";
const double PTPlacement::UNASSIGNED_POSTQ = EGriceLab::nan;
const double PTPlacement::UNASSIGNED_DIST = EGriceLab::nan;
const double PTPlacement::UNASSIGNED_RATIO = EGriceLab::nan;

void PTPlacement::calcQValues(vector<PTPlacement>& places, PRIOR_TYPE type) {
	if(places.empty())
		return;

	/* explore all placements */
	VectorXd ppPlace(places.size()); /* posterior logP at placement */
	map<string, double> ppTaxon; /* posterior logP at taxon */
	double ppTaxNorm = infV; /* log(0) */

	VectorXd::Index i = 0;
	for(vector<PTPlacement>::const_iterator placement = places.begin(); placement != places.end(); ++placement) {
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
	for(vector<PTPlacement>::size_type i = 0; i < places.size(); ++i) {
		double q = EGriceLab::Math::p2q(1 - p(i));
		places[i].qPlace = q > MAX_Q ? MAX_Q : q;
	}

	/* calculate qTaxonomy */
	for(vector<PTPlacement>::iterator placement = places.begin(); placement != places.end(); ++placement) {
		double q = EGriceLab::Math::p2q(1 - ::exp(ppTaxon[placement->getTaxonName()] - ppTaxNorm));
		placement->qTaxon = q > MAX_Q ? MAX_Q : q;
	}
}

/**
 * calculate prior probability at log-scale
 * @param place  a placement
 * @param type  prior type
 * @param h  base height of this placement (for cNode)
 * @return  log prior always no greater than 0
 */
double PTPlacement::logPriorPr(PRIOR_TYPE type) const {
	double logP;
	switch(type) {
	case UNIFORM:
		logP = -0;
		break;
	case HEIGHT:
		logP = -(annoDist - wnr + height);
		break;
	}
	return logP;
}

string alignSeq(const BandedHMMP7& hmm, const CSFMIndex& csfm, const PrimarySeq& read,
		int seedLen, int seedRegion, BandedHMMP7::align_mode mode,
		int& csStart, int& csEnd) {
	const DegenAlphabet* abc = hmm.getNuclAbc();
	const int K = hmm.getProfileSize();
	const int L = read.length();

	BandedHMMP7::ViterbiScores seqVscore(K, L); // construct an empty reusable score
	vector<BandedHMMP7::ViterbiAlignPath> seqVpaths; // construct an empty list of VPaths
	BandedHMMP7::ViterbiAlignTrace seqVtrace; // construct an empty VTrace

	int regionLen = seedRegion < read.length() ? seedRegion : read.length(); /* search region */
	/* find seed in 5' */
	for(int seedFrom = 0; seedFrom < regionLen - seedLen + 1; ++seedFrom) {
		PrimarySeq seed(abc, read.getId(), read.subseq(seedFrom, seedLen));
		const CSLoc& loc = csfm.locateOne(seed.getSeq());
		if(loc.start > 0 && loc.end > 0) /* a read seed located */ {
//					fprintf(stderr, "start:%d end:%d from:%d to:%d  CSLen:%d CS:%s\n", loc.start, loc.end, seedFrom + 1, seedFrom + seedLen, loc.CS.length(), loc.CS.c_str());
			seqVpaths.push_back(hmm.buildAlignPath(loc, seedFrom + 1, seedFrom + seedLen)); /* seed_from and seed_to are 1-based */
			break; /* only one 5'-seed necessary */
		}
	}
	/* find seed in 3', if requested */
	if(mode == BandedHMMP7::GLOBAL) {
		for(int seedTo = read.length() - 1; seedTo - seedLen + 1 >= read.length() - regionLen; --seedTo) {
			PrimarySeq seed(abc, read.getId(), read.subseq(seedTo - seedLen + 1, seedLen));
			const CSLoc& loc = csfm.locateOne(seed.getSeq());
			if(loc.start > 0 && loc.end > 0) /* a read seed located */ {
//						fprintf(stderr, "start:%d end:%d from:%d to:%d  CSLen:%d CS:%s\n", loc.start, loc.end, seedTo - seedLen + 2, seedTo + 1, loc.CS.length(), loc.CS.c_str());
				seqVpaths.push_back(hmm.buildAlignPath(loc, seedTo - seedLen + 2, seedTo + 1)); /* seed_from and seed_to are 1-based */
				break; /* only one 3'-seed necessary */
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

	/* find seqStart and seqEnd */
	csStart = hmm.getCSLoc(seqVtrace.alnStart);
	csEnd = hmm.getCSLoc(seqVtrace.alnEnd);

	return hmm.buildGlobalAlign(read, seqVscore, seqVtrace);
}

vector<PTLoc> getSeed(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, double maxDist) {
	vector<PTLoc> locs; /* candidate locations */
	/* get potential placement locations based on pDist to observed or inferred sequences */
	for(vector<PTUnrooted::PTUNodePtr>::size_type i = 0; i < ptu.numNodes(); ++i) {
		PTUnrooted::PTUNodePtr node = ptu.getNode(i);
		if(node->isRoot())
			continue;
		double pDist = SeqUtils::pDist(node->getSeq(), seq, start, end);
		if(pDist <= maxDist)
			locs.push_back(PTLoc(node, pDist));
	}
	return locs;
}

vector<PTPlacement> estimateSeq(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, const vector<PTLoc>& locs) {
	vector<PTPlacement> places;
	for(vector<PTLoc>::const_iterator loc = locs.begin(); loc < locs.end(); ++loc) {
		const PTUnrooted::PTUNodePtr& cNode = loc->node;
		const PTUnrooted::PTUNodePtr& pNode = cNode->getParent();
		double cDist = loc->dist;
		double pDist = SeqUtils::pDist(pNode->getSeq(), seq, start, end);
		double ratio = cDist / (cDist + pDist);
		if(::isnan(ratio)) // unable to estimate the ratio
			ratio = 0.5;
//		cerr << "Estimating at " << cNode->getId() << " cDist: " << cDist << " pDist: " << pDist << " ratio: " << ratio << endl;
		/* estimate the placement */
		double wnr;
		double loglik = ptu.estimateSeq(seq, cNode, pNode, start, end, ratio, wnr);
		places.push_back(PTPlacement(cNode, pNode, ratio, wnr, loglik));
	}
	return places;
}

/** Get accurate placement for a seq given the estimated placements */
vector<PTPlacement>& placeSeq(const PTUnrooted& ptu, const DigitalSeq& seq, int start, int end,
		vector<PTPlacement>& places) {
	/* accurate placement using estimated values */
	for(vector<PTPlacement>::iterator place = places.begin(); place != places.end(); ++place) {
		double ratio0 = place->ratio;
		double wnr0 = place->wnr;
		double loglik0 = place->loglik;

//		cerr << "Estimated placement ratio0: " << ratio0 << " wnr: " << wnr0 << " loglik: " << loglik0 << endl;
		PTUnrooted subtree = ptu.copySubTree(place->cNode, place->pNode);
		const PTUnrooted::PTUNodePtr& v = subtree.getNode(0);
		const PTUnrooted::PTUNodePtr& u = subtree.getNode(1);
		double w0 = subtree.getBranchLength(u, v);

		double loglik = subtree.placeSeq(seq, u, v, start, end, ratio0, wnr0);
		const PTUnrooted::PTUNodePtr& r = subtree.getNode(2);
		const PTUnrooted::PTUNodePtr& n = subtree.getNode(3);

		/* update placement info */
		double wnr = subtree.getBranchLength(n, r);
		double wur = subtree.getBranchLength(u, r);
		double wvr = w0 - wur;
		place->ratio = wur / w0;
//		cerr << "delta loglik: " << (loglik - loglik0) << endl;
		place->height = ptu.getHeight(place->cNode) + wur;
		place->annoDist = wvr <= wur ? wvr + wnr : wur + wnr;
		/* update other placement info */
		place->wnr = wnr;
		place->loglik = loglik;
	} /* end each candidate placement */

	return places;
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

} /* End namespace EGriceLab */
