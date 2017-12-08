#!/usr/bin/env bash

# Input
inLoc=$1
oLoc=$2
factor=$3

# Creating MEME fasta input
outLoc=$oLoc$factor"_TEMP/"
mkdir -p $outLoc
gunzip -c $inLoc"ChIPseq."$factor".conservative.train.narrowPeak.gz" > $outLoc"a.bed"
python /work/eg474423/eg474423_Projects/trunk/Dream/exp/meme_chip/narrowPeakToCenter.py "50" $outLoc"a.bed" $outLoc"b.bed"
sort -k1,1 -k2,2n $outLoc"b.bed" > $outLoc"c.bed"
fastaFromBed -fi "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa" -bed $outLoc"c.bed" -fo $outLoc"a.fa"

# Applying MEME-chip
memeLoc=$oLoc$factor"/"
mkdir -p $memeLoc
meme-chip -oc $memeLoc -db "/work/eg474423/eg474423_Projects/trunk/Dream/exp/meme_chip/JASPAR_CORE_2016_vertebrates.meme" -nmeme "600" -meme-mod "zoops" -meme-minw "5" -meme-maxw "15" -meme-nmotifs "5" -meme-p "4" $outLoc"a.fa"
rm -rf $outLoc

# Fetching motifs
# TODO


