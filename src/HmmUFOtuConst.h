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
 * HmmUFOtuConst.h
 *
 *  Created on: Jun 3, 2015
 *      Author: zhengqi
 */
#include<string>
#include<iostream>
#include<limits>
#include<cassert>

#ifndef HMMUFOTUCONST_H_
#define HMMUFOTUCONST_H_

namespace EGriceLab {
using std::string;
using std::istream;
using std::ostream;

/* includes for some program-wide constants and definitions */
const std::string progName = "HmmUFOtu";
const std::string progVersion = "v1.02";

const double inf = std::numeric_limits<double>::infinity();
const double infV = -inf;
const double nan = std::numeric_limits<double>::quiet_NaN();

const string MSA_FILE_SUFFIX = ".msa";
const string CSFM_FILE_SUFFIX = ".csfm";
const string HMM_FILE_SUFFIX = ".hmm";
const string SUB_MODEL_FILE_SUFFIX = ".sm";
const string PHYLOTREE_FILE_SUFFIX = ".ptu";

const string GZIP_FILE_SUFFIX = ".gz";
const string BZIP2_FILE_SUFFIX = ".bz2";

/* stand-alone functions */
/**
 * Compare two versions
 */
int cmpVersion(const string& ver1, const string& ver2);

/**
 * Read progName from an input stream
 */
istream& readProgName(istream& in, string& name);

/**
 * Read progVersion from an input stream
 */
istream& readProgVersion(istream& in, string& version);

/**
 * Write progName to an output stream
 */
ostream& writeProgName(ostream& out, const string& name);

/**
 * Write progVersion to an output stream
 */
ostream& writeProgVersion(ostream& out, const string& version);

/**
 * show program and package version
 */
void printVersion(const string& prog);

} /* end namespace EGriceLab */
#endif /* HMMUFOTUCONST_H_ */
