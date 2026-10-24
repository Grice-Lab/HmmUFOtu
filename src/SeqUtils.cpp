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
 * SeqUtils.cpp
 *
 *  Created on: May 10, 2017
 *      Author: zhengqi
 */

#include <cassert>
#include <algorithm>
#include "StringUtils.h"
#include "SeqUtils.h"

namespace EGriceLab {
namespace HmmUFOtu {

const string SeqUtils::FASTA_FMT = "fasta";
const string SeqUtils::FASTQ_FMT = "fastq";

const char *SeqUtils::FASTA_FILE_EXTENSIONS[] = { "fasta", "fas", "fa", "fna" };
const char *SeqUtils::FASTQ_FILE_EXTENSIONS[] = { "fastq", "fq" };


double SeqUtils::pDist(const DigitalSeq& seq1, const DigitalSeq& seq2,
		DigitalSeq::size_type start, DigitalSeq::size_type end) {
	assert(seq1.getAbc() == seq2.getAbc());
	assert(seq1.length() == seq2.length());

	DigitalSeq::size_type d = 0;
	DigitalSeq::size_type N = 0;
	const DegenAlphabet* abc = seq1.getAbc();
	for(DigitalSeq::size_type i = start; i <= end; ++i) {
		DigitalSeq::value_type b1 = seq1[i];
		DigitalSeq::value_type b2 = seq2[i];
		if(abc->isSymbol(b1) && abc->isSymbol(b2)) { // is a symbol
			N++;
			if(b1 != b2)
				d++;
		}
	}
	return static_cast<double>(d) / N;
}

double SeqUtils::pDist(const string& seq1, const string& seq2,
		string::size_type start, string::size_type end) {
	assert(seq1.length() == seq2.length());
	string::size_type d = 0;
	for(string::size_type i = start; i <= end; ++i)
		if(seq1[i] != seq2[i])
			d++;
	return static_cast<double>(d) / (end - start + 1);
}

double SeqUtils::pDist(const string& seq1, const string& seq2,
		const DegenAlphabet* abc, string::size_type start,
		string::size_type end) {
	assert(seq1.length() == seq2.length());
	string::size_type d = 0;
	string::size_type N = 0;
	for(string::size_type i = start; i <= end; ++i) {
		char c1 = seq1[i];
		char c2 = seq2[i];
		if(abc->isSymbol(c1) && abc->isSymbol(c2)) { /* only count non-gaps */
			N++;
			if(!abc->isMatch(c1, c2))
				d++;
		}
	}
	return static_cast<double>(d) / N;
}

double SeqUtils::pDist(const string& seq1, const DigitalSeq& seq2, size_t start,
		size_t end) {
	assert(seq1.length() == seq2.length());
	const DegenAlphabet* abc = seq2.getAbc();
	size_t d = 0;
	size_t N = 0;
	for(size_t i = start; i <= end; ++i) {
		char c = seq1[i];
		int8_t b = seq2[i];
		if(abc->isSymbol(c) && abc->isSymbol(b)) { /* ignore gaps */
			N++;
			if(!abc->isMatch(c, b))
				d++;
		}
	}
	return static_cast<double>(d) / N;
}

bool SeqUtils::isFastaFileExt(const string& fn) {
	for(const char **ext = FASTA_FILE_EXTENSIONS;
			ext != FASTA_FILE_EXTENSIONS + sizeof(FASTA_FILE_EXTENSIONS) / sizeof(*FASTA_FILE_EXTENSIONS); ++ext)
		if(StringUtils::endsWith(fn, *ext))
			return true;
	return false;
}

bool SeqUtils::isFastqFileExt(const string& fn) {
	for(const char **ext = FASTQ_FILE_EXTENSIONS;
			ext != FASTQ_FILE_EXTENSIONS + sizeof(FASTQ_FILE_EXTENSIONS) / sizeof(*FASTQ_FILE_EXTENSIONS); ++ext)
		if(StringUtils::endsWith(fn, *ext))
			return true;
	return false;
}

string SeqUtils::guessSeqFileFormat(const string& fn) {
	if(isFastaFileExt(fn))
		return FASTA_FMT;
	else if(isFastqFileExt(fn))
		return FASTQ_FMT;
	else
		return "";
}

} /* namespace HmmUFOtu */
} /* namespace EGriceLab */
