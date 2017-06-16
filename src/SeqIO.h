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
 *
 *  Created on: Jul 23, 2015
 *      Author: zhengqi
 */

#ifndef SEQIO_H_
#define SEQIO_H_

#include <fstream>
#include "PrimarySeq.h"

namespace EGriceLab {

using std::string;
using std::ifstream;
using std::ofstream;
/**
 * A class to handle IO operation for PrimarySeq of various format and
 */
class SeqIO {
public:
	/* nested class and enums */
	enum Mode { READ, WRITE };

	/* constructors */
	/** default constructor */
	SeqIO() : abc(NULL) {  }

	/**
	 * Construct a SeqIO object with given info
	 */
	SeqIO(const string& filename, const string& alphabet, const string& format, Mode mode = READ, int maxLine = DEFAULT_MAX_LINE);

	/**
	 * Construct a SeqIO object with given info
	 */
	SeqIO(const string& filename, const DegenAlphabet* abc, const string& format, Mode mode = READ, int maxLine = DEFAULT_MAX_LINE);

	/* Getters and Setters */
	const string& getFilename() const {
		return filename;
	}

	const string& getFormat() const {
		return format;
	}

	const Mode getMode() const {
		return mode;
	}

	int getMaxLine() const {
		return maxLine;
	}

	void setMaxLine(int maxLine) {
		this->maxLine = maxLine;
	}

	/* member methods */
	/** open a new SeqIO, close old one if necessary */
	void open(const string& filename, const DegenAlphabet* abc, const string& format, Mode mode = READ, int maxLine = DEFAULT_MAX_LINE);

	/** open a new SeqIO, close old one if necessary */
	void open(const string& filename, const string& alphabet, const string& format, Mode mode = READ, int maxLine = DEFAULT_MAX_LINE) {
		open(filename, AlphabetFactory::getAlphabetByName(alphabet), format, mode, maxLine);
	}

	/**
	 * close underlying iostreams explicitly
	 */
	void close() {
		if(in.is_open())
			in.close();
		if(out.is_open())
			out.close();
	}

	/** test whether this SeqIO is open */
	bool is_open() const {
		return mode == READ ? in.is_open() : out.is_open();
	}

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
	string filename;
	string format;
	const DegenAlphabet* abc;
	Mode mode;
	int maxLine;

	ifstream in;
	ofstream out;
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
	if(!in.is_open())
		throw std::ios_base::failure("Unable to read from an unopened file");
	if(format == "fasta")
		return nextFastaSeq();
	else if(format == "fastq")
		return nextFastqSeq();
	else
		return PrimarySeq(abc, "", "");
}

inline void SeqIO::writeSeq(const PrimarySeq& seq) {
	if(!out.is_open())
		throw std::ios_base::failure("Unable to write to an unopened file");
	if(format == "fasta")
		writeFastaSeq(seq);
	else if(format == "fastq")
		writeFastqSeq(seq);
	else { } /* do nothing */
}

} /* namespace EGriceLab */

#endif /* SEQIO_H_ */
