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
 * HmmUFOtu_main.h
 *  Type definitions for hmmufotu core algorithms
 *  Created on: Jul 10, 2017
 *  	Since: v1.1
 *      Author: zhengqi
 */

#ifndef SRC_HMMUFOTU_MAIN_H_
#define SRC_HMMUFOTU_MAIN_H_

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include "HmmUFOtu.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::string;
using std::vector;

/**
 * A JSON Placement type for holding an HmmUFOtu placement result
 */
struct JPlace {
	/* constructors */
	/** default constructor */
	JPlace() {  }

	/** construct a JPlacement from PTPlacment info */
	JPlace(int edgeID, string readName, double edgeLen, double ratio,
			double loglik, double annoDist, double q);

	/* member fields */
	int edgeID;
	string readName;
	double likelihood;
	double like_ratio;
	double distal_length;
	double proximal_length;
	double pendant_length;

	/* static member fields */
	static const int MAX_Q = 250; /* maximum allowed Q value */
};

/** Align seq using banded HMM algorithm, returns an HmmAlignment */
BandedHMMP7::HmmAlignment alignSeq(const BandedHMMP7& hmm, const CSFMIndex& csfm, const PrimarySeq& read,
		int seedLen, int seedRegion, BandedHMMP7::align_mode mode);

/** Align seq using traditional HMM algorithm, returns an HmmAlignment */
BandedHMMP7::HmmAlignment alignSeq(const BandedHMMP7& hmm, const PrimarySeq& read);

/**
 * Get seed placement locations by checking p-dist between a given seq and observed/inferred seq of nodes
 * @param ptu  PTUnrooted tree to be used
 * @param seq  sequence to be placed
 * @param start  0-based start
 * @param end  0-based end
 * @param maxDiff  maximum allowed p-Distance difference
 * @param maxHeight  maximum allowed height of nodes to place
 * @return  a vector of PTPlacement sorted by the p-dist
 */
vector<PTUnrooted::PTLoc> getSeed(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, double maxDiff, double maxHeight);

/** Get estimated placement for a seq at given locations */
vector<PTUnrooted::PTPlacement> estimateSeq(const PTUnrooted& ptu, const DigitalSeq& seq,
		const vector<PTUnrooted::PTLoc>& locs, const string& method);

/**
 * filter estimated placement by removing bad placement with estimated loglik lower than the best placement
 * @param places  a vector of placements
 * @param maxError  maximum error of log-liklihood allowed compared to the best placement
 * @return  the modified vector of placements sorted by their loglike decreasingly
 */
vector<PTUnrooted::PTPlacement>& filterPlacements(vector<PTUnrooted::PTPlacement>& places, double maxError);

/** Get accurate placement for a seq given the estimated placements */
vector<PTUnrooted::PTPlacement>& placeSeq(const PTUnrooted& ptu, const DigitalSeq& seq,
		vector<PTUnrooted::PTPlacement>& places, double maxHeight = inf);

/** calculate Q-values using a given prior type */
void calcQValues(vector<PTUnrooted::PTPlacement>& places, PTUnrooted::PRIOR_TYPE type);

/** get alignment identity, as fraction of non-gap characters in the alignment part */
double alignIdentity(const DegenAlphabet* abc, const string& align, int start, int end);

/** get profile-HMM identity, as fraction of non-gap characters in HMM profile sites */
double hmmIdentity(const BandedHMMP7& hmm, const string& align, int start, int end);

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_HMMUFOTU_MAIN_H_ */
