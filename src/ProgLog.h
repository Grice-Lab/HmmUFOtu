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
using std::cerr;
using std::cout;
using std::endl;

class ProgLog {
public:
	/* nested enum */
	enum LogLevel {
		LOG_NOTHING,
//		LOG_CRITICAL,
		LOG_ERROR,
		LOG_WARNING,
		LOG_INFO,
		LOG_DEBUG
	};

	/* constructors */
	/**
	 * construct a ProgLog using given output stream
	 */
	explicit ProgLog(ostream& out, LogLevel level = LOG_NOTHING)
	: out(out), level(level) {  }
	/* disable copy and assignment constructors */
private:
	ProgLog(const ProgLog& other);
	ProgLog& operator=(const ProgLog& other);

public:
	/* member methods */
	LogLevel getLevel() const {
		return level;
	}

	void setLevel(LogLevel level) {
		this->level = level;
	}

	template <typename T>
	ProgLog& operator<<(const T& value) {
		if(VERBOSE_LEVEL >= level)
			out << value;
		return *this;
	}

	ProgLog& operator<<(ostream& (*pf)(ostream&)) {
		if(VERBOSE_LEVEL >= level)
			out << *pf;
		return *this;
	}

private:
	ostream& out; /* underlying output stream */
	LogLevel level;


};

/* namespace static variables */
extern ProgLog errorLog;
extern ProgLog warningLog;
extern ProgLog infoLog;

} /* namespace EGriceLab */

#endif /* SRC_PROGLOG_H_ */
