/*
 * PTLoc.h
 *
 *  Created on: Feb 8, 2018
 *      Author: zhengqi
 */

#ifndef SRC_PTLOC_H_
#define SRC_PTLOC_H_

#include "PhyloTreeUnrooted.h"

namespace EGriceLab {
namespace HmmUFOtu {


/**
 * A "seed" Phylogenetic tree loc to store potential placement location information
 */
struct PTLoc {
	/* constructors */
	/** default constructor */
	PTLoc() {  }

	/** construct a location using given node and distance */
	PTLoc(int start, int end, const PTUnrooted::PTUNodePtr& node, double dist)
	: start(start), end(end), node(node), dist(dist)
	{  }

	/* non-member functions */
	friend bool operator<(const PTLoc& lhs, const PTLoc& rhs);

	/* member fields */
	int start;
	int end;
	PTUnrooted::PTUNodePtr node;
	double dist;
};


inline bool operator<(const PTLoc& lhs, const PTLoc& rhs) {
	return lhs.dist < rhs.dist;
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */

#endif /* SRC_PTLOC_H_ */
