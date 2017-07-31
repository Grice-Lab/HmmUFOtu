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
 * SeqIO.h
 *  A class for popular sequence file IO
 *  Created on: Jul 23, 2015
 *      Author: zhengqi
 */

#ifndef SEQIO_H_
#define SEQIO_H_

#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include "PrimarySeq.h"

namespace EGriceLab {

using std::string;
using std::istream;
using std::ostream;
using std::streambuf;
using std::ifstream;
using std::ofstream;
/**
 * A class to handle IO operation for PrimarySeq of various format and
 */
class SeqIO {
public:
	/* constructors */
	/** default constructor, do nothing */
	SeqIO() : abc(NULL), in(NULL), out(NULL) {  }

	/**
	 * Construct a SeqIO object in READ mode with given info
	 */
	SeqIO(istream* in, const DegenAlphabet* abc, const string& format, int maxLine = DEFAULT_MAX_LINE);

	/**
	 * Construct a SeqIO object in WRITE mode with given info
	 */
	SeqIO(ostream* out, const DegenAlphabet* abc, const string& format, int maxLine = DEFAULT_MAX_LINE);

	/** destructor, do nothing */
	virtual ~SeqIO() {  }

public:

	/* Getters and Setters */
	const string& getFormat() const {
		return format;
	}

	int getMaxLine() const {
		return maxLine;
	}

	void setMaxLine(int maxLine) {
		this->maxLine = maxLine;
	}

	/* member methods */
	/** set the input to a given a new istream, will not close the old one */
	void reset(istream* in, const DegenAlphabet* abc, const string& format, int maxLine = DEFAULT_MAX_LINE);

	/** set the out to a given a new ostream, will not close the old one */
	void reset(ostream* out, const DegenAlphabet* abc, const string& format, int maxLine = DEFAULT_MAX_LINE);

	/**
	 * test whether this file has next PrimarySeq
	 * @return true if everything is good and has symbol indicating nextSeq exists
	 */
	bool hasNext();

	/**
	 * Get next PrimarySeq, if possible
	 * @return PrimarySeq, if hasNext is true, otherwise return an empty seq with everything empty
	 * @throw std::ios_base::failure if nextSeq not available or other IO exception
	 */
	PrimarySeq nextSeq();

	/**
	 * Write a seq to the output
	 * @param seq  a PrimarySeq
	 * @throw std::ios_base::failure if any IO exception
	 */
	void writeSeq(const PrimarySeq& seq);

private:
	/* Disable copy and assign constructors */
	SeqIO(const SeqIO& other);
	SeqIO& operator=(const SeqIO& other);

	/**
	 * Get next PrimarySeq in fasta format, if possible
	 * @return PrimarySeq, if hasNext is true, otherwise return an empty seq with everything empty
	 * @throw std::ios_base::failure if nextSeq not available or other IO exception
	 */
	PrimarySeq nextFastaSeq();

	/**
	 * Get next PrimarySeq in fasta format, if possible
	 * @return PrimarySeq, if hasNext is true, otherwise return an empty seq with everything empty
	 * @throw std::ios_base::failure if nextSeq not available or other IO exception
	 */
	PrimarySeq nextFastqSeq();

	/**
	 * test whether this file has next PrimarySeq in fasta format
	 * @return true if everything is good and has symbol indicating nextSeq exists
	 */
	bool hasNextFasta();

	/**
	 * test whether this file has next PrimarySeq in fastq format
	 * @return true if everything is good and has symbol indicating nextSeq exists
	 */
	bool hasNextFastq();

	/**
	 * Write a seq to the output in fasta format
	 * @param seq  a PrimarySeq
	 * @throw std::ios_base::failure if any IO exception
	 */
	void writeFastaSeq(const PrimarySeq& seq);

	/**
	 * Write a seq to the output in fastq format,
	 * with maxLine restricted
	 * @param seq  a PrimarySeq
	 * @param maxLine  max characters in a line, set to -1 for limits
	 * @throw std::ios_base::failure if any IO exception
	 */
	void writeFastqSeq(const PrimarySeq& seq);

private:
	/** member fields */
	string format;
	const DegenAlphabet* abc;
	int maxLine;

	istream* in; /* input */
	ostream* out; /* output */

	/* static members */
	static const char fastaHead = '>';
	static const char fastqHead = '@';
	static const int DEFAULT_MAX_LINE = 60;
	static const char fastqSep = '+';
};

inline bool SeqIO::hasNext() {
	if(format == "fasta")
		return hasNextFasta();
	else if(format == "fastq")
		return hasNextFastq();
	return false;
}

inline PrimarySeq SeqIO::nextSeq() {
	if(format == "fasta")
		return nextFastaSeq();
	else if(format == "fastq")
		return nextFastqSeq();
	else
		return PrimarySeq(abc, "", "");
}

inline void SeqIO::writeSeq(const PrimarySeq& seq) {
	if(format == "fasta")
		writeFastaSeq(seq);
	else if(format == "fastq")
		writeFastqSeq(seq);
	else { } /* do nothing */
}

} /* namespace EGriceLab */

#endif /* SEQIO_H_ */
