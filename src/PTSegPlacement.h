/*
 * PTSegPlacement.h
 *
 *  Created on: Feb 9, 2018
 *      Author: zhengqi
 */

#ifndef SRC_PTSEGPLACEMENT_H_
#define SRC_PTSEGPLACEMENT_H_

#include "PTPlacement.h"

namespace EGriceLab {
namespace HmmUFOtu {

struct PTSegPlacement: public PTPlacement {
public:
	/* constructors */
	/** default constructor */
	/** default constructor */
	PTSegPlacement() : start(0), end(0), segStart(0), segEnd(0), ratio(nan), wnr(nan), loglik(nan),
	height(nan), annoDist(nan), qPlace(nan), qTaxon(nan)
	{  }

	/** construct a SegPlacement with basic info and optionally auxilary info */
	PTSegPlacement(int start, int end, int segStart, int segEnd,
			const PTUnrooted::PTUNodePtr& cNode, const PTUnrooted::PTUNodePtr& pNode,
			double ratio, double wnr, double loglik,
			double height = 0, double annoDist = 0,
			double qPlace = 0, double qTaxonomy = 0)
	: PTPlacement(start, end, cNode, pNode, ratio, wnr, loglik, height, annoDist, qPlace, qTaxonomy),
	  segStart(segStart), segEnd(segEnd)
	{  }

	virtual ~PTSegPlacement() {  }

	/* member methods */
	bool isPlaced() const {
		return treeLoglik.rows() > 0;
	}

	double segLoglik(int segStart, int segEnd) const {
		return treeLoglik.segment(segStart, segEnd - segStart + 1).sum();
	}

	/* member fields */
	int segStart;
	int segEnd;
	VectorXd treeLoglik;
};

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_PTSEGPLACEMENT_H_ */
