#!/bin/bash

INPUT="70_otus"
DBNAME="gg_70_otus"
SMTYPE="GTR"
DB="${DBNAME}_${SMTYPE}_dG"
SRCPATH="../src"
DATAPATH="../data"

echo $DB

echo "Constructing a test database ..."
$SRCPATH/hmmufotu-build ${INPUT}.fasta ${INPUT}.tree -a ${INPUT}_taxonomy.txt -n $DB -s $SMTYPE -V -v
if [ $? == 0 ]
	then
		echo "$DB constructed successfully"
	else
		echo "Failed to constructing $DB"
		rm -f ${DB}.*
		exit 1
fi

echo "Testing MSA IO ..."
./MSAIO_test ${DB}.msa ${DB}.2.msa
if [ $? == 0 ]
	then
		echo "MSA IO passed"
	else
		echo "MSA IO failed"
		rm -f ${DB}.*
		exit 1
fi

echo "Testing CSFM IO ..."
./FMIO_test ${DB}.csfm ${DB}.2.csfm
if [ $? == 0 ]
	then
		echo "CSFM IO passed"
	else
		echo "CSFM IO failed"
		rm -f ${DB}.*
		exit 1
fi

echo "Testing bHMM IO ..."
./bHmm_IO_test ${DB}.hmm ${DB}.2.hmm
if [ $? == 0 ]
	then
		echo "bHMM IO passed"
	else
		echo "bHMM IO failed"
		rm -f ${DB}.*
		exit 1
fi

echo "Testing PTU IO ..."
./PTU_IO_test ${DB}.ptu ${DB}.2.ptu
if [ $? == 0 ]
	then
		echo "PTU IO passed"
	else
		echo "PTU IO failed"
		rm -f ${DB}.*
		exit 1
fi

rm -f ${DB}.*
