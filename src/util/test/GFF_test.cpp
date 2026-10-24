/*
 * GFF_test.cpp
 *
 *  Created on: Aug 7, 2018
 *      Author: zhengqi
 */

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include "GFF.h"
using namespace std;
using namespace EGriceLab::UCSC;

int main() {
	GFF gffGene("1", "transcribed_unprocessed_pseudogene", "gene",
			11869, 14409, GFF::INVALID_SCORE, '+', GFF::INVALID_FRAME,
			GFF::attr_map { {"gene_id", "ENSG00000223972"},
							{"gene_name",  "DDX11L1"},
							{"gene_source", "havana"},
							{"gene_biotype", "transcribed_unprocessed_pseudogene"} });

	GFF gffRna(GFF::GFF3, "ctg123", ".", "mRNA",
			1300, 9000, GFF::INVALID_SCORE, '+', GFF::INVALID_FRAME,
			"ID=mrna0001;Name=sonichedgehog");

	ostringstream gtfout, gff3out;

	/* formatted IO tests */
	gffGene.write(gtfout, GFF::GTF);
	if(gtfout.bad()) {
		cerr << "Unable to write gffGene in GTF format: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	string gtfoutstr = gtfout.str();
	cout << "gtfout content:" << endl << gtfoutstr << endl;

	gff3out << gffRna;
	if(gff3out.bad()) {
		cerr << "Unable to write gffRna in GFF3 format: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}
	string gff3outstr = gff3out.str();
	cout << "gff3out content:" << endl << gff3outstr << endl;

	GFF gffGeneN, gffRnaN;
	istringstream gtfin(gtfoutstr);
	istringstream gff3in(gff3outstr);
	gffGeneN.read(gtfin, GFF::GTF);
	if(gtfin.bad()) {
		cerr << "Unable to read gffGene in GTF format: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	gff3in >> gffRnaN;
	if(gff3in.bad()) {
		cerr << "Unable to read gffRna in GFF3 format: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(gffGeneN != gffGene) {
		cerr << "Read gffGene doesn't equal original copy" << endl;
		return EXIT_FAILURE;
	}

	if(gffRnaN != gffRna) {
		cerr << "Read gffRna doesn't equal original copy" << endl;
		return EXIT_FAILURE;
	}

	/* binary IO tests */
	ostringstream out;
	gffGene.save(out);
	if(out.bad()) {
		cerr << "Unable to save gffGene: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	gffRna.save(out);
	if(out.bad()) {
		cerr << "Unable to save gffRna: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	gffGeneN.clear();
	gffRnaN.clear();
	istringstream in(out.str());
	gffGeneN.load(in);
	if(in.bad()) {
		cerr << "Unable to load gffGene: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	gffRnaN.load(in);
	if(in.bad()) {
		cerr << "Unable to load gffRna: " << ::strerror(errno) << endl;
		return EXIT_FAILURE;
	}

	if(gffGeneN != gffGene) {
		cerr << "Loaded gffGene doesn't equal original copy" << endl;
		return EXIT_FAILURE;
	}

	if(gffRnaN != gffRna) {
		cerr << "Loaded gffRna doesn't equal original copy" << endl;
		return EXIT_FAILURE;
	}

	return 0;
}
