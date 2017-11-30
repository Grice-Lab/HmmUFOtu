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
 * TSVScanner.h
 *  A Scanner to parse TSVRecord from an istream
 *  Created on: Jul 12, 2017
 *      Author: zhengqi
 */

#ifndef SRC_TSVSCANNER_H_
#define SRC_TSVSCANNER_H_

#include <map>
#include <iostream>
#include <cstdarg>
#include "TSVRecord.h"

using std::string;
using std::map;
using std::ifstream;
using std::ofstream;

namespace EGriceLab {

class TSVScanner {
public:
	/* constructors */
	/** construct a TSVScanner with given input and aux info */
	TSVScanner(istream& in, bool hasHeader = false,
			const string& sep = DEFAULT_SEP, char quote = DEFAULT_QUOTE);

	/** destructor, do nothing */
	virtual ~TSVScanner() {  }

	/* disable copy and assign operators */
private:
	TSVScanner(const TSVScanner& other);
	TSVScanner& operator=(const TSVScanner& other);

public:
	/** member methods */
	/** test whether this TSVScanner has a header */
	bool hasHeader() const {
		return header != NULL && !header->empty();
	}

	/**
	 * test whether this file has next record
	 * @return true if everything is good and has additional lines
	 */
	bool hasNext();

	/**
	 * Get next TSVRecord, if possible
	 * @return PrimarySeq, if hasNext is true, otherwise return an empty seq with everything empty
	 * @throw std::ios_base::failure if nextRecord not available or other IO exception
	 */
	TSVRecord nextRecord();

private:
	istream& in;
	TSVRecord::TSVHeaderPtr header;
	string sep;
	char quote;

	/** static fields */
public:
	static const char COMMENT_CHAR = '#';
	static const string DEFAULT_SEP;
	static const char DEFAULT_QUOTE = '\0';

	/** internal methods */
private:
	void parseHeader();
};


} /* namespace EGriceLab */

#endif /* SRC_TSVSCANNER_H_ */
