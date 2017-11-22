/*
 * TSVIO.cpp
 *
 *  Created on: Jul 12, 2017
 *      Author: zhengqi
 */

#include "TSVScanner.h"

#include <cstdio> /* EOF def */
#include <cassert>

namespace EGriceLab {
using namespace std;

const string TSVScanner::DEFAULT_SEP = "\t";

TSVScanner::TSVScanner(istream& in, bool hasHeader, const string& sep, char quote) :
		in(in), sep(sep), quote(quote) {
	if(hasHeader)
		parseHeader();
}

bool TSVScanner::hasNext() {
	return in.peek() != EOF;
}

TSVRecord TSVScanner::nextRecord() {
	string line;
	std::getline(in, line);
	return TSVRecord(line, header, sep, quote);
}

void TSVScanner::parseHeader() {
	string line;
	while(!hasHeader() && in.peek() != EOF) {
		std::getline(in, line);
		if(line[0] == COMMENT_CHAR) // a comment line near the header
			continue;
		/* construct a new header */
		header.reset(new TSVRecord::TSVHeader(line, sep)); /* construct a new header from header line */
	}
}

} /* namespace EGriceLab */
