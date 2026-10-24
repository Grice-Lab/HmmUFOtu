/*
 * SAMtools_test.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <cstdio>
#include "BAMheader.h"
#include "HTSindex.h"
#include "BAM.h"
#include "SAMfile.h"

using namespace std;
using namespace EGriceLab::SAMtools;

int main(int argc, char *argv[]) {
	/* construct a BAMheader */
	BAMheader::targetMap chrMap {
		{ "chr1", 1000000 },
		{ "chr2", 500000  }
	};

	BAMheader header(chrMap);
	cerr << "BAMheader constructed" << endl;
	header.setHDTag("SQ", "unsorted");

	header.addTag("@PG", argv[0]);
	cerr << "BAMheader tag added" << endl;

	/* construct BAM records */
	const uint32_t readLen = 100;

	BAM bam1("r1", 0, header.getIndex("chr1"), 100000, 0, BAM::encodeCigar("100M"), string(readLen, 'A'), BAM::qual_str(readLen, 20));
	BAM bam2("r2", 0, header.getIndex("chr1"), 100000, 0, BAM::encodeCigar("50M1X49M"), string(readLen, 'C'), BAM::qual_str(readLen, 20));
	BAM bam3("r3", 0, header.getIndex("chr1"), 100000, 0, BAM::encodeCigar("10S90M"), string(readLen, 'G'), BAM::qual_str(readLen, 20));
	BAM bam4("r4", 0, header.getIndex("chr1"), 100000, 0, BAM::encodeCigar("40=2X8=1X49="), string(readLen, 'T'), BAM::qual_str(readLen, 20));
	BAM bamNull("r5", string(readLen, 'N'), BAM::qual_str(readLen, 2));
	cerr << "bam constructed" << endl;

	bam1.setAux("XI", 10);
	bam2.setAux("XD", 10.0);
	bam3.setAux("XZ", "str");
	int32_t arr[4] { 0, 1, 2, 3 };
	bam4.setAux("XA", 'i', 4, arr);
	cerr << "bam aux added" << endl;

	/* construct SAM and BAM files */
	string samFn = string(argv[0]) + ".sam";
	string bamFn = string(argv[0]) + ".bam";
	SAMfile samOut(samFn, "w");
	SAMfile bamOut(bamFn, "wb");
	cerr << "SAMfile constructed" << endl;

	samOut.setHeader(header);
	bamOut.setHeader(header);
	cerr << "SAMfile header set" << endl;

	/* IO tests */
	samOut.writeHeader();
	if(samOut.write(bam1) == -1) {
		cerr << "Failed to write bam1 to '" << samFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(samOut.write(bam2) == -1) {
		cerr << "Failed to write bam2 to '" << samFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(samOut.write(bam3) == -1) {
		cerr << "Failed to write bam3 to '" << samFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(samOut.write(bam4) == -1) {
		cerr << "Failed to write bam4 to '" << samFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(samOut.write(bamNull) == -1) {
		cerr << "Failed to write bamNull to '" << samFn << "'" << endl;
		return EXIT_FAILURE;
	}
	cerr << "SAM IO tested" << endl;

	bamOut.writeHeader();
	if(bamOut.write(bam1) == -1) {
		cerr << "Failed to write bam1 to '" << bamFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(bamOut.write(bam2) == -1) {
		cerr << "Failed to write bam2 to '" << bamFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(bamOut.write(bam3) == -1) {
		cerr << "Failed to write bam3 to '" << bamFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(bamOut.write(bam4) == -1) {
		cerr << "Failed to write bam4 to '" << bamFn << "'" << endl;
		return EXIT_FAILURE;
	}
	if(bamOut.write(bamNull) == -1) {
		cerr << "Failed to write bamNull to '" << bamFn << "'" << endl;
		return EXIT_FAILURE;
	}
	cerr << "BAM IO tested" << endl;

	samOut.close();
	bamOut.close();
	int status = SAMfile::buildIndex(bamFn);
	if(status != 0) {
		cerr << "Failed to build index for '" << bamFn << "': " << status << endl;
		return EXIT_FAILURE;
	}
	/* remove temporary files */
//	remove(samFn.c_str());
//	remove(bamFn.c_str());
//	remove((bamFn + HTSindex::DEFAULT_INDEX_SUFFIX).c_str());
}
