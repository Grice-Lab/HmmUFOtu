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
	/* constructors */
	/**
	 * construct a ProgLog using given output stream
	 */
	explicit ProgLog(ostream& out, LOG_LEVEL level = LOG_NOTHING)
	: out(out), level(level) {  }
	/* disable copy and assignment constructors */
private:
	ProgLog(const ProgLog& other);
	ProgLog& operator=(const ProgLog& other);

public:
	/* member methods */
	LOG_LEVEL getLevel() const {
		return level;
	}

	void setLevel(LOG_LEVEL level) {
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
	LOG_LEVEL level;
};

/* namespace static variables */
extern ProgLog errorLog;
extern ProgLog warningLog;
extern ProgLog infoLog;
extern ProgLog debugLog;

} /* namespace EGriceLab */

#endif /* SRC_PROGLOG_H_ */
