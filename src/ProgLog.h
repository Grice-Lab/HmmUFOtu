/*
 * ProgLog.h
 *  A customized template program-logging class that use formatted output using the '<<' operator
 *  Created on: Dec 14, 2016
 *      Author: zhengqi
 */

#ifndef SRC_PROGLOG_H_
#define SRC_PROGLOG_H_

#include <iostream>
#include "HmmUFOtuEnv.h"

namespace EGriceLab {
using std::ostream;
using std::streambuf;
using std::cerr;
using std::cout;
using std::endl;

void UPDATE_LOGS();

inline void DISABLE_ALL() {
	VERBOSE_LEVEL = LOG_NOTHING;
	UPDATE_LOGS();
}

inline void INCREASE_LEVEL(int increment = 1) {
	VERBOSE_LEVEL += increment;
	UPDATE_LOGS();
}

inline void DECREASE_LEVEL(int decrement = 1) {
	VERBOSE_LEVEL -= decrement;
	UPDATE_LOGS();
}

/* namespace static variables */
extern ostream errorLog;
extern ostream warningLog;
extern ostream infoLog;
extern ostream debugLog;

} /* namespace EGriceLab */

#endif /* SRC_PROGLOG_H_ */
