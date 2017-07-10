#!/bin/zsh

# Export
export CC=gcc # GCC
export JAVA_TOOL_OPTIONS=-Xmx12000m
export FPICFLAGS=-fPIC
module load gcc &> /dev/null
module load python &> /dev/null

export PATH=$PATH:/hpcwork/izkf/bin/
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/hpcwork/izkf/lib64/python2.7/site-packages

# Input
bedFileName=$1
pfmFileNameList=$2
genomeFileName=$3
footprintBamFileName=$4
dnaseBamFileName=$5
outputFileName=$6

# Execution
python ./fullMpbsFootprintTc_Final.py $bedFileName $pfmFileNameList $genomeFileName $footprintBamFileName $dnaseBamFileName $outputFileName

