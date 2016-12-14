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
extern int VERBOSE_LEVEL; /* DEFAULT VERBOSE LEVEL */

enum LOG_LEVEL {
	LOG_NOTHING,
	LOG_ERROR,
	LOG_WARNING,
	LOG_INFO,
	LOG_DEBUG
};

void DISABLE_VERBOSE_MESSAGES();

void ENABLE_ERROR();

void ENABLE_WARNING();

void ENABLE_INFO();

}

#endif /* SRC_HMMUFOTUENV_H_ */
