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

void DISABLE_VERBOSE_MESSAGES() {
	VERBOSE_LEVEL = LOG_NOTHING;
	errorLog.rdbuf(NULL);
	warningLog.rdbuf(NULL);
	infoLog.rdbuf(NULL);
	debugLog.rdbuf(NULL);
}

void ENABLE_ERROR() {
	if(VERBOSE_LEVEL < LOG_ERROR) {
		VERBOSE_LEVEL = LOG_ERROR;
		errorLog.rdbuf(std::cerr.rdbuf());
	}
}

void ENABLE_WARNING() {
	if(VERBOSE_LEVEL < LOG_WARNING) {
		VERBOSE_LEVEL = LOG_WARNING;
		errorLog.rdbuf(std::cerr.rdbuf());
		warningLog.rdbuf(std::cerr.rdbuf());
	}
}

void ENABLE_INFO() {
	if(VERBOSE_LEVEL < LOG_INFO) {
		VERBOSE_LEVEL = LOG_INFO;
		errorLog.rdbuf(std::cerr.rdbuf());
		warningLog.rdbuf(std::cerr.rdbuf());
		infoLog.rdbuf(std::cerr.rdbuf());
	}
}

void ENABLE_DEBUG() {
	if(VERBOSE_LEVEL < LOG_DEBUG) {
		VERBOSE_LEVEL = LOG_DEBUG;
		errorLog.rdbuf(std::cerr.rdbuf());
		warningLog.rdbuf(std::cerr.rdbuf());
		infoLog.rdbuf(std::cerr.rdbuf());
		debugLog.rdbuf(std::cerr.rdbuf());
	}
}
} /* namespace EGriceLab */
