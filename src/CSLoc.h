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
	CSLoc() : start(0), end(0) {  }

	/**
	 * Construct a CSLoc at given loc
	 */
	CSLoc(int start, int end, const string& CS = "")
		: start(start), end(end), CS(CS)
	{  }

	/** member methods */
	bool isValid() const {
		return start > 0 && start < end && CS.length() > end - start;
	}

	bool isValid(int from, int to) const {
		return isValid() && CS.length() > to - from;
	}

	int start; // CS start
	int end;   // CS end
	string CS; // CS string
};

} /* namespace EGriceLab */

#endif /* SRC_CSLOC_H_ */
