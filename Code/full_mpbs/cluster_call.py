
# Import
import os
import sys

bedFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/test_regions.blacklistfiltered.bed"
bedFileName = "/work/ig440396/project/Dream/exp/full_mpbs/chrom.sizes.hg19.bed"
genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
motifListFileName = "/work/ig440396/project/Dream/exp/motifs/pfm_hocomoco/info.txt"
outLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs/hocomoco_test/"

motifListFile = open(motifListFileName,"r")

# Cell Loop
for line in motifListFile:

  ll = line.strip().split("\t")
  mLoc = "/".join(motifListFileName.split("/")[:-1])+"/"
  motif = ll[0]
  if(ll[1] == "."): continue
  pfmFileNameList = ",".join([mLoc+e+".pwm" for e in ll[1].split(",")])
  outputFileName = outLoc+motif+".bed"

  # Execution cluster
  myL = "BS_"+motif
  clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
  clusterCommand += "-W 120:00 -M 12000 -S 100 -P izkf -R \"select[hpcwork]\" ./pipeline.zsh "
  clusterCommand += bedFileName+" "+pfmFileNameList+" "+genomeFileName+" "+outputFileName
  #os.system(clusterCommand)
  print(clusterCommand)

motifListFile.close()


