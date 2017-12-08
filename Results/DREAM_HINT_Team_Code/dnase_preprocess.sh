#!/usr/bin/env bash

# Input
originalLoc=$1
finalBamLoc=$2
cell=$3

# Merge BAM
cd $originalLoc
samtools merge $cell"_DNase_merged.bam" "./DNASE."$cell"."*".bam"
samtools sort $cell"_DNase_merged.bam" -o $cell"_DNase.bam"
rm $cell"_DNase_merged.bam"
samtools index $cell"_DNase.bam"

# Quality Assessment
samtools rmdup -sS $cell"_DNase.bam" $cell"_remdup.bam"
samtools view -bq 20 $cell"_remdup.bam" > $cell"_DNase.bam"
samtools index $cell"_DNase.bam"
rm $cell"_remdup.bam"

# Move BAM to final location
mkdir -p $finalBamLoc$cell
cd $finalBamLoc$cell
mv $originalLoc$cell"_DNase.bam" "."
mv $originalLoc$cell"_DNase.bam.bai" "."

# Apply MACS
mkdir -p "Peaks"
macs2 callpeak -t $cell"_DNase.bam" -n "Peaks" --outdir "./Peaks/" -f BAM -g "hs" --nomodel --nolambda --keep-dup auto --call-summits

cd ..


