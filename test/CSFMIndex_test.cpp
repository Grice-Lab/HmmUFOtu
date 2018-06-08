/*
 * CSFMIndex_test.cpp
 *
 *  Created on: Jun 8, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include "HmmUFOtu_common.h"
#include "HmmUFOtu_hmm.h"

using namespace std;
using namespace EGriceLab::HmmUFOtu;

int main() {
	string alnSeq =
			string(">seq1\nATCA-ctg\n") +
			string(">seq2\nATCCGG-T\n") +
			string(">seq3\nATCGC-GT\n") +
			string(">seq4\nATCTCGG-\n");

	istringstream in(alnSeq);

	MSA msa;
	msa.loadMSA(in, "fasta");
	if(in.bad()) {
		cerr << "Failed to load MSA: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	CSFMIndex csfm;
	csfm.build(msa);
	cout << "CSFM-Index build" << endl;

	string pat = "ATC";
	/* test count */
	int32_t count = csfm.count(pat);
	cout << "Found " << count << " matches of " << pat << " in MSA: " << alnSeq;
	if(count != 4)
		return EXIT_FAILURE;
	/* test locate */
	vector<CSLoc> locs = csfm.locate(pat);
	for(vector<CSLoc>::const_iterator loc = locs.begin(); loc != locs.end(); ++loc) {
		cout << "Found matched CSLoc: " << loc->start << "-" << loc->end << endl;
		if(!(loc->start == 1 && loc->end == 3))
			return EXIT_FAILURE;
	}
	/* test locateFirst */
	CSLoc loc = csfm.locateFirst(pat);
	cout << "Found first matched CSLoc: " << loc.start << "-" << loc.end << endl;
	if(!(loc.start == 1 && loc.end == 3))
		return EXIT_FAILURE;

	/* test locateRandom */
	loc = csfm.locateFirst(pat);
	cout << "Found random matched CSLoc: " << loc.start << "-" << loc.end << endl;
	if(!(loc.start == 1 && loc.end == 3))
		return EXIT_FAILURE;
}
