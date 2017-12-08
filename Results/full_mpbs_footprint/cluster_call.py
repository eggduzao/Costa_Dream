
# Import
import os
import sys
from glob import glob


#bedFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/test_regions.blacklistfiltered.bed"
#bedFileName = "/work/ig440396/project/Dream/exp/full_mpbs/chrom.sizes.hg19.bed"
footprintFiles = "/hpcwork/izkf/projects/dream_tfbs/exp/footprints/results/"
genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
motifListFileName = "/work/ig440396/project/Dream/exp/motifs/pfm_hocomoco/info_fpr.txt"
outLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_tc/hocomoco/"

inFileList = glob(footprintFiles+"*.bed")

tfCellTableFileName = "../tf_cell_table.txt"
tfCellTableFile = open(tfCellTableFileName,"r")
list={}
for l in tfCellTableFile:
  l=l.strip("\n")
  l=l.split("\t")
  tf=l[0]
  for cell in l[2].split("."):
    list[(tf,cell)]=""
  for cell in l[3].split("."):
    list[(tf,cell)]=""


for inF in inFileList:
  cell = inF.split("/")[-1].split(".")[0]

  motifListFile = open(motifListFileName,"r")
  #motif loop
  for line in motifListFile:

    try:
      ll = line.strip().split("\t")
      mLoc = "/".join(motifListFileName.split("/")[:-1])+"/"
      motif = ll[0]
      x=list[(motif,cell)]
    
      bitscore = ll[2]
      if(ll[1] == "."): continue
      pfmFileNameList = ",".join([mLoc+e+".pwm" for e in ll[1].split(",")])
      outputFileName = outLoc+motif+"."+cell+".bed"

      # Execution cluster
      myL = "BS_"+motif
      clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
      clusterCommand += "-W 120:00 -M 12000 -S 100 -P izkf -R \"select[hpcwork]\" ./pipeline.zsh "
      clusterCommand += inF+" "+pfmFileNameList+" "+genomeFileName+" "+outputFileName+" "+bitscore
      os.system(clusterCommand)
      #print(clusterCommand)
    except:
      continue

  motifListFile.close()


