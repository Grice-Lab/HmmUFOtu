#!/bin/bash

# basic info
DB="gg_70_otus_GTR"
SRCPATH="../src"
DATAPATH="../data"

# simulating info
SIMFILE="${DB}_sim.fasta"
SIMNUM=100
SIMSEED=0

# assigning info
ASSIGNFILE="${DB}_sim_assign.txt"
CHIMERAFILE="${DB}_sim_chimera.txt"

# OTU info
OTUFILE="${DB}_sim_OTU.txt"
OTULIST="${DB}_sim_OTU.list"
OTUALIGN="${DB}_sim_OTU.fasta"
OTUTREE="${DB}_sim_OTU.tree"

# OTU manipulation
SUBSETFILE="${DB}_sim_OTU_subset.txt"
SUBSETN=20
NORMFILE="${DB}_sim_OTU_norm.txt"
MERGEDFILE="${DB}_sim_OTU_merged.txt"
MERGEDTREE="${DB}_sim_OTU_merged.tree"

echo "Using database $DB"

echo "Generating simulated reads ..."
$SRCPATH/hmmufotu-sim $DB $SIMFILE -N $SIMNUM -S $SIMSEED -v
if [ $? == 0 ]
	then
		echo "simulated reads generated"
	else
		echo "Failed to generate simulated reads"
		exit 1
fi

echo "Running taxonomy assignment ..."
$SRCPATH/hmmufotu $DB $SIMFILE -o $ASSIGNFILE -v
if [ $? == 0 ]
	then
		echo "taxonomy assignment file generated"
	else
		echo "Failed to generate assignment file"
		exit 1
fi 

echo "Running taxonomy assignment with chimera checking enabled ..."
$SRCPATH/hmmufotu $DB $SIMFILE -o $ASSIGNFILE -v -C --chimera-out $CHIMERAFILE
if [ $? == 0 ]
then
echo "taxonomy assignment file generated"
else
echo "Failed to generate assignment file"
exit 1
fi


echo "Summarizing OTU table ..."
$SRCPATH/hmmufotu-sum $DB $ASSIGNFILE -o $OTUFILE -r $OTULIST -c $OTUALIGN -t $OTUTREE -v
if [ $? == 0 ]
	then
		echo "OTU table generated"
	else
		echo "Failed to generate OTU table and related files"
		exit 1
fi 

echo "Subsetting OTU table ..."
$SRCPATH/hmmufotu-subset $OTUFILE $SUBSETFILE -s $SUBSETN -v
if [ $? == 0 ]
	then
		echo "OTU table subset"
	else
		echo "Failed to subset OTU table"
		exit 1
fi 

echo "Normalizing OTU table ..."
$SRCPATH/hmmufotu-norm $OTUFILE $NORMFILE -v
if [ $? == 0 ]
	then
		echo "OTU table normalized"
	else
		echo "Failed to normalize OTU table"
		exit 1
fi 

echo "Merging OTU table ..."
$SRCPATH/hmmufotu-merge $SUBSETFILE $NORMFILE -o $MERGEDFILE --db $DB -t $MERGEDTREE -v
if [ $? == 0 ]
	then
		echo "OTU tables merged"
	else
		echo "Failed to merge OTU tables"
		exit 1
fi 

rm -f ${DB}*

