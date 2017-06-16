/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
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
