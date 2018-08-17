/*
 * GFF_test.cpp
 *
 *  Created on: Aug 7, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "GFF.h"
using namespace std;
using namespace EGriceLab::UCSC;
using namespace EGriceLab;

int main() {
	const string gtfinstr = "1\ttranscribed_unprocessed_pseudogene\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972\"; gene_name \"DDX11L1\"; gene_source \"havana\"; gene_biotype \"transcribed_unprocessed_pseudogene\"";
	const string gff3instr = "ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mrna0001;Name=sonichedgehog";

	cout << "gtfin:" << endl << gtfinstr << endl;
	cout << "gff3in:" << endl << gff3instr << endl;

	istringstream gtfin(gtfinstr);
	istringstream gff3in(gff3instr);

	ostringstream gtfout, gff3out;

	GFF gtf(GFF::GTF);

	gtfin >> gtf;
	gtfout << gtf;
	if(gtfout.str() != gtfinstr) {
		cerr << "Unmatched gtfout:" << endl << gtfout.str() << endl;
		return 1;
	}

	GFF gff3(GFF::GFF3);
	gff3in >> gff3;
	gff3out << gff3;
	if(gff3out.str() != gff3instr) {
		cerr << "Unmatched gff3out:" << endl << gff3out.str() << endl;
		return 1;
	}

	return 0;
}
