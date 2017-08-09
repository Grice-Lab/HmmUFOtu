/*
 * CSLoc.h
 *
 *  Created on: Aug 9, 2017
 *      Author: zhengqi
 *      Since: v1.1
 */

#ifndef SRC_CSLOC_H_
#define SRC_CSLOC_H_

#include <string>

namespace EGriceLab {
using std::string;
/**
 * A public class for describing a region on the consensus seq (CS)
 */
struct CSLoc {
	/* constructors */
	/**
	 * Default constructor, do nothing
	 */
	CSLoc() {  }

	/**
	 * Construct a CSLoc at given loc
	 */
	CSLoc(int start, int end, const string& CS = "")
		: start(start), end(end), CS(CS)
	{  }

	int start; // CS start
	int end;   // CS end
	string CS; // CS string
};

} /* namespace EGriceLab */

#endif /* SRC_CSLOC_H_ */
