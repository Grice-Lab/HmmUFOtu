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
