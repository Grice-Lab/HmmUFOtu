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

# jplace info
OTUFILE="${DB}_sim_assign.jplace"

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

echo "Converting assignment into jplace format ..."
$SRCPATH/hmmufotu-jplace $DB $ASSIGNFILE -o $OTUFILE -sm -V -a -v
if [ $? == 0 ]
	then
		echo "jplace file converted"
	else
		echo "Failed to convert assignments into jplace file"
		exit 1
fi 
