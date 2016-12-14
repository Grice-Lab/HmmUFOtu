/*
 * HmmUFOtuEnv.cpp
 *
 *  Created on: Jul 15, 2016
 *      Author: zhengqi
 */

#include "HmmUFOtuEnv.h"

namespace EGriceLab {

int VERBOSE_LEVEL = 2; /* DEFAULT VERBOSE LEVEL */

void DISABLE_VERBOSE_MESSAGES() {
	VERBOSE_LEVEL = LOG_NOTHING;
}

void ENABLE_ERROR() {
	if(VERBOSE_LEVEL < LOG_ERROR)
		VERBOSE_LEVEL = LOG_ERROR;
}

void ENABLE_WARNING() {
	if(VERBOSE_LEVEL < LOG_WARNING)
		VERBOSE_LEVEL = LOG_WARNING;
}

void ENABLE_INFO() {
	if(VERBOSE_LEVEL < LOG_INFO)
		VERBOSE_LEVEL = LOG_INFO;
}

}
