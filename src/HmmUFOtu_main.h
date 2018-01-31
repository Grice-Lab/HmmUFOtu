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
#include <iostream>
#include <algorithm>
#include <cmath>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include "HmmUFOtu.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::string;
using std::vector;
using std::ostream;

/** prior probability method */
enum PRIOR_TYPE {
	UNIFORM,
	HEIGHT
};

/** A "seed" Phylogenetic tree loc to store potential placement location information */
struct PTLoc {
	/* constructors */
	/** construct a location using given node and distance */
	PTLoc(const PTUnrooted::PTUNodePtr& node, double dist)
	: node(node), dist(dist)
	{  }

	/* non-member functions */
	friend bool operator<(const PTLoc& lhs, const PTLoc& rhs);

	PTUnrooted::PTUNodePtr node;
	double dist;
};

/** A candidate Phylogenetic Tree Placement to store (partially) placement information */
struct PTPlacement {
//	/** default constructor */
	PTPlacement() : ratio(nan), wnr(nan), loglik(nan),
			height(nan), annoDist(nan), qPlace(nan), qTaxon(nan)
	{ }

	/** construct a placement with basic info and optionally auxilary info */
	PTPlacement(const PTUnrooted::PTUNodePtr& cNode, const PTUnrooted::PTUNodePtr& pNode,
			double ratio, double wnr, double loglik,
			double height = 0, double annoDist = 0,
			double qPlace = 0, double qTaxonomy = 0)
	: cNode(cNode), pNode(pNode), ratio(ratio), wnr(wnr), loglik(loglik),
	  height(height), annoDist(annoDist), qPlace(qPlace), qTaxon(qTaxonomy)
	{ }

	/** member methods */
	long getTaxonId() const {
		if(cNode != NULL && pNode != NULL)
			return ratio <= 0.5 ? cNode->getId() : pNode->getId();
		else
			return UNASSIGNED_TAXONID;
	}

	string getTaxonName() const {
		if(cNode != NULL && pNode != NULL)
			return ratio <= 0.5 ? cNode->getAnno() : pNode->getAnno();
		else
			return UNASSIGNED_TAXONNAME;
	}

	string getId() const {
		if(cNode != NULL && pNode != NULL)
			return boost::lexical_cast<string> (cNode->getId()) + "->" + boost::lexical_cast<string> (cNode->getParent()->getId());
		else
			return UNASSIGNED_ID;
	}

	/** calculate prior probability of a placement given a prior type in log-scale */
	double logPriorPr(PRIOR_TYPE type) const;

	/** calculate prior proability of a placement given a prior type */
	double priorPr(PRIOR_TYPE type) const {
		return ::exp(logPriorPr(type));
	}

	/** static member methods */
	/** calculate Q-values using a given prior type */
	static void calcQValues(vector<PTPlacement>& places, PRIOR_TYPE type);

	/** non-member functions */
	friend bool compareByLoglik(const PTPlacement& lhs, const PTPlacement& rhs);
	friend bool compareByQTaxon(const PTPlacement& lhs, const PTPlacement& rhs);
	friend bool compareByQPlace(const PTPlacement& lhs, const PTPlacement& rhs);
	friend ostream& operator<<(ostream& out, const PTPlacement& place);

	/** member fields */
	PTUnrooted::PTUNodePtr cNode;
	PTUnrooted::PTUNodePtr pNode;
	double ratio; /* placement ratio */
	double wnr;   /* new branch length */
	double loglik;
	double annoDist;
	double height;
	double qPlace;
	double qTaxon;

	/** static member fields */
	static const int MAX_Q = 250; /* maximum allowed Q value */
	static const long UNASSIGNED_TAXONID = -1;
	static const string UNASSIGNED_TAXONNAME;
	static const double UNASSIGNED_LOGLIK;
	static const string UNASSIGNED_ID;
	static const double UNASSIGNED_POSTQ;
	static const double UNASSIGNED_DIST;
	static const double UNASSIGNED_RATIO;
};

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

/** Align seq with hmm and csfm, returns the alignment and update csStart and csEnd */
string alignSeq(const BandedHMMP7& hmm, const CSFMIndex& csfm, const PrimarySeq& read,
		int seedLen, int seedRegion, BandedHMMP7::align_mode mode,
		int& csStart, int& csEnd);

/**
 * Get seed placement locations by checking p-dist between a given seq and observed/inferred seq of nodes
 * @param ptu  PTUnrooted
 * @param seq  sequence to be placed
 * @param start  0-based start
 * @param end  0-based end
 * @param maxDiff  maximum allowed p-Distance difference
 * @return  a vector of PTLoc sorted by the p-dist
 * */
vector<PTLoc> getSeed(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, double maxDiff);

/** Get estimated placement for a seq at given locations */
vector<PTPlacement> estimateSeq(const PTUnrooted& ptu, const DigitalSeq& seq,
		int start, int end, const vector<PTLoc>& locs, const string& method);

/** Get accurate placement for a seq given the estimated placements */
vector<PTPlacement>& placeSeq(const PTUnrooted& ptu, const DigitalSeq& seq, int start, int end,
		vector<PTPlacement>& places);

/** get alignment identity, as fraction of non-gap characters in the alignment part */
double alignIdentity(const DegenAlphabet* abc, const string& align, int start, int end);

/** get profile-HMM identity, as fraction of non-gap characters in HMM profile sites */
double hmmIdentity(const BandedHMMP7& hmm, const string& align, int start, int end);

inline bool operator<(const PTLoc& lhs, const PTLoc& rhs) {
	return lhs.dist < rhs.dist;
}

inline ostream& operator<<(ostream& out, const PTPlacement& place) {
	out << place.getId() << "\t" << place.ratio << "\t"
			<< place.getTaxonId() << "\t" << place.getTaxonName() << "\t"
			<< place.annoDist << "\t" << place.loglik << "\t"
			<< place.qPlace << "\t" << place.qTaxon;
	return out;
}

inline bool compareByLoglik(const PTPlacement& lhs, const PTPlacement& rhs) {
	return lhs.loglik < rhs.loglik;
}

inline bool compareByQPlace(const PTPlacement& lhs, const PTPlacement& rhs) {
	return lhs.qPlace < rhs.qPlace;
}

inline bool compareByQTaxon(const PTPlacement& lhs, const PTPlacement& rhs) {
	return lhs.qTaxon < rhs.qTaxon;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_HMMUFOTU_MAIN_H_ */
