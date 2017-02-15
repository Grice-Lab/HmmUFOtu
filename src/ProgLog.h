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

void DISABLE_VERBOSE_MESSAGES();

void ENABLE_ERROR();

void ENABLE_WARNING();

void ENABLE_INFO();

void ENABLE_DEBUG();

/* namespace static variables */
extern ostream errorLog;
extern ostream warningLog;
extern ostream infoLog;
extern ostream debugLog;

} /* namespace EGriceLab */

#endif /* SRC_PROGLOG_H_ */
