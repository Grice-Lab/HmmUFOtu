/*
 * BandedHMMCommons.h
 *	Contains basic public class for BandedHMM algorithm
 *  Created on: Aug 3, 2015
 *      Author: zhengqi
 */

#ifndef BANDEDHMMCOMMONS_H_
#define BANDEDHMMCOMMONS_H_
#include <string>

namespace EGriceLab {
using std::string;
/**
 * A public class for describing a region on the concensus seq (CS)
 */
struct CSLoc {
	/* constructors */
	/**
	 * Construct a CSLoc at given loc
	 */
	CSLoc(int start, int end, const string& CS, int rank = 1)
		: start(start), end(end), CS(CS), rank(rank) { }

	int start;
	int end;
	string CS;
	int rank;
};

}

#endif /* BANDEDHMMCOMMONS_H_ */
