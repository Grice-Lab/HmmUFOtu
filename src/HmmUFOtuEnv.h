/*
 * HmmUFOtuEnv.h
 *  Per-application environment variables
 *  Created on: Nov 18, 2015
 *      Author: zhengqi
 */

#ifndef SRC_HMMUFOTUENV_H_
#define SRC_HMMUFOTUENV_H_
#include <string>
#include <iostream>

namespace EGriceLab {
/* per-application variables */
static int verbose = 0;

void printVerboseInfo(const std::string& message) {
	if(verbose > 0)
		std::cerr << message << std::endl;
}

}

#endif /* SRC_HMMUFOTUENV_H_ */
