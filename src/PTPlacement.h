/*
 * PTPlacement.h
 *
 *  Created on: Feb 8, 2018
 *      Author: zhengqi
 */

#ifndef SRC_PTPLACEMENT_H_
#define SRC_PTPLACEMENT_H_

#include <iostream>
#include <boost/lexical_cast.hpp>
#include "PhyloTreeUnrooted.h"

namespace EGriceLab {
namespace HmmUFOtu {

using std::ostream;

/**
 * A candidate Phylogenetic Tree Placement to store (partially) placement information
 */
struct PTPlacement {
	/* nested types */
	/** prior probability enum */
	enum PRIOR_TYPE {
		UNIFORM,
		HEIGHT
	};

	/* constructors */
//	/** default constructor */
	PTPlacement() : start(0), end(0), ratio(nan), wnr(nan), loglik(nan),
			height(nan), annoDist(nan), qPlace(nan), qTaxon(nan)
	{ }

	/** construct a placement with basic info and optionally auxilary info */
	PTPlacement(int start, int end,
			const PTUnrooted::PTUNodePtr& cNode, const PTUnrooted::PTUNodePtr& pNode,
			double ratio, double wnr, double loglik,
			double height = 0, double annoDist = 0,
			double qPlace = 0, double qTaxonomy = 0)
	: start(start), end(end), cNode(cNode), pNode(pNode),
	  ratio(ratio), wnr(wnr), loglik(loglik),
	  height(height), annoDist(annoDist), qPlace(qPlace), qTaxon(qTaxonomy)
	{  }

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
			return boost::lexical_cast<string> (cNode->getId()) + "->" + boost::lexical_cast<string> (pNode->getId());
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
	int start;
	int end;
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

	static const string TSV_HEADER;
};

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

#endif /* SRC_PTPLACEMENT_H_ */
