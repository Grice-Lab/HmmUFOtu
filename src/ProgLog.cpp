/*
 * ProgLog.cpp
 *
 *  Created on: Dec 14, 2016
 *      Author: zhengqi
 */

#include "ProgLog.h"

namespace EGriceLab {

using namespace std;

ostream errorLog(VERBOSE_LEVEL >= LOG_ERROR ? std::cerr.rdbuf() : NULL);
ostream warningLog(VERBOSE_LEVEL >= LOG_WARNING ? std::cerr.rdbuf() : NULL);
ostream infoLog(VERBOSE_LEVEL >= LOG_INFO ? std::cerr.rdbuf() : NULL);
ostream debugLog(VERBOSE_LEVEL >= LOG_DEBUG ? std::cerr.rdbuf() : NULL);

void UPDATE_LOGS() {
	errorLog.rdbuf(VERBOSE_LEVEL >= LOG_ERROR ? std::cerr.rdbuf() : NULL);
	warningLog.rdbuf(VERBOSE_LEVEL >= LOG_WARNING ? std::cerr.rdbuf() : NULL);
	infoLog.rdbuf(VERBOSE_LEVEL >= LOG_INFO ? std::cerr.rdbuf() : NULL);
	debugLog.rdbuf(VERBOSE_LEVEL >= LOG_DEBUG ? std::cerr.rdbuf() : NULL);
}

} /* namespace EGriceLab */
