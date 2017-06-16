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
 * SeqIO.cpp
 *
 *  Created on: Jul 23, 2015
 *      Author: zhengqi
 */
#include <fstream>
#include <cctype>
#include "SeqIO.h"

namespace EGriceLab {
using namespace std;

SeqIO::SeqIO(const string& filename, const string& alphabet, const string& format, Mode mode, int maxLine) :
	filename(filename), abc(AlphabetFactory::getAlphabetByName(alphabet)), format(format), mode(mode), maxLine(maxLine) {
	/* check format support */
	if(!(format == "fasta" || format == "fastq"))
		throw invalid_argument("Unsupported file format '" + format + "'");
	/* register exceptions */
//	in.exceptions(std::ifstream::badbit);
//	out.exceptions(std::ofstream::badbit);

	/* open files */
	if(mode == READ)
		in.open(filename.c_str());
	else
		out.open(filename.c_str());
}

SeqIO::SeqIO(const string& filename, const DegenAlphabet* abc, const string& format, Mode mode, int maxLine) :
	filename(filename), abc(abc), format(format), mode(mode), maxLine(maxLine) {
	/* check format support */
	if(!(format == "fasta" || format == "fastq"))
		throw invalid_argument("Unsupported file format '" + format + "'");
	/* register exceptions */
//	in.exceptions(std::ifstream::badbit);
//	out.exceptions(std::ofstream::badbit);

	/* open files */
	if(mode == READ)
		in.open(filename.c_str());
	else
		out.open(filename.c_str());
}

void SeqIO::open(const string& filename, const DegenAlphabet* abc, const string& format, Mode mode, int maxLine) {
	/* close old resources */
	close();
	/* check format support */
	if(!(format == "fasta" || format == "fastq"))
		throw invalid_argument("Unsupported file format '" + format + "'");
	/* replace values */
	this->filename = filename;
	this->abc = abc;
	this->format = format;
	this->mode = mode;
	this->maxLine = maxLine;

	/* open files */
	if(mode == READ)
		in.open(filename.c_str());
	else
		out.open(filename.c_str());
}

bool SeqIO::hasNextFasta() {
	char c = in.peek();
	return c != EOF && c == fastaHead;
}

bool SeqIO::hasNextFastq() {
	char c = in.peek();
	return c != EOF && c == fastqHead;
}

PrimarySeq SeqIO::nextFastaSeq() {
	string id, seq, desc;
	char tag;
	string line;
	tag = in.get();
	if(tag != fastaHead)
		throw ios_base::failure("inputfile " + filename + " is not a valid FASTA format");

	in >> id; // read the next word as id
	while(::isspace(in.peek()) && in.peek() != '\n') // ignore non-newline white spaces
		in.get();
	getline(in, desc); // read the remaining as desc, if any
	while(in.peek() != EOF && in.peek() != fastaHead) {
		getline(in, line);
		seq += line;
	}
	return PrimarySeq(abc, id, seq, desc);
}

PrimarySeq SeqIO::nextFastqSeq() {
	string id, seq, desc, qual;
	char tag;
	string line;
	tag = in.get();
	if(tag != fastqHead)
		throw ios_base::failure("inputfile " + filename + " is not a valid FASTQ format");
	in >> id; // read the next word as id
	while(::isspace(in.peek()) && in.peek() != '\n') // ignore non-newline white spaces
		in.get();
	getline(in, desc); // read the remaining as description
	getline(in, seq);  // read seq line
	getline(in, line); // ignore sep line
	getline(in, qual); // read qual line
	return PrimarySeq(abc, id, seq, desc, qual);
}

void SeqIO::writeFastaSeq(const PrimarySeq& seq) {
	out << fastaHead << seq.getId() << (!seq.getDesc().empty() ? " " + seq.getDesc() : "") << endl;
	if(maxLine > 0) {
		const char* seqPtr = seq.getSeq().c_str();
		for(size_t i = 0, r = seq.length(); i < seq.length(); i += maxLine, r -= maxLine) {
			out.write(seqPtr + i, r >= maxLine ? maxLine : r); /* use unformated write for performance */
			out.put('\n'); // do not flush for faster performance
		}
	}
	else
		out << seq.getSeq() << endl;
}

void SeqIO::writeFastqSeq(const PrimarySeq& seq) {
	out << fastqHead << seq.getId() << (!seq.getDesc().empty() ? " " + seq.getDesc() : "") << endl;
	out << seq.getSeq() << endl;
	out << fastqSep << endl << seq.getQual() << endl;
}

} /* namespace EGriceLab */

