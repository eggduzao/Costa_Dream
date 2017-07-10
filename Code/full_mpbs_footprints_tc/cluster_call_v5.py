
# Import
import os
import sys
from glob import glob

# Iterating on the challenge level
challengeList = ["train", "test", "ladder"]
for challenge in challengeList:

  # Challenge Parameters
  if(challenge == "ladder"):
    challengeLabel = "L"
    challengeCol = 2
  elif(challenge == "train"):
    challengeLabel = "T"
    challengeCol = 1
  else:
    challengeLabel = "F"
    challengeCol = 3

  # Iterating on TF-CELL table
  tfCellTableFileName = "/work/eg474423/Costa_Dream/exp/tf_cell_table.txt"
  tfCellTableFile = open(tfCellTableFileName,"r")
  tfCellTableFile.readline()
  for line in tfCellTableFile:

    # Line parsing
    ll = line.strip().split("\t")
    if(ll[challengeCol] == "."): continue
    factor = ll[0]
    if((challenge == "train") and (factor != "CTCF" and factor != "E2F1" and factor != "EGR1" and factor != "FOXA1" and factor != "FOXA2" and factor != "GABPA" and factor != "HNF4A" and factor != "JUND" and factor != "MAX" and factor != "NANOG" and factor != "REST" and factor != "TAF1")): continue
    if((challenge == "test") and (factor != "CTCF" and factor != "E2F1" and factor != "EGR1" and factor != "FOXA1" and factor != "FOXA2" and factor != "GABPA" and factor != "HNF4A" and factor != "JUND" and factor != "MAX" and factor != "NANOG" and factor != "REST" and factor != "TAF1")): continue
    if((challenge == "ladder") and (factor != "CTCF" and factor != "EGR1")): continue
    trainCellList = ll[1].split(",")
    testCellList = ll[challengeCol].split(",")

    # Creating PFM list based on training cells
    pfmLoc = "/work/eg474423/Costa_Dream/exp/motifs/pfm_hocomoco_meme/"
    pfmInfoFileName = pfmLoc+"info.txt"
    pfmInfoFile = open(pfmInfoFileName,"r")
    pfmList = []
    for line2 in pfmInfoFile:
      ll2 = line2.strip().split("\t")
      if(ll2[0] == factor):
        pfmList = [pfmLoc+e+".pwm" for e in ll2[1].split(",")]
        break
    pfmInfoFile.close()

    # Iterating on test cells
    for cell in testCellList:

      # Iterating on chromosomes
      bedLoc = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/split_"+challenge+"_chrom/"

      chr1L = ["chr1."+str(e) for e in range(1,11)]
      chr2L = ["chr2."+str(e) for e in range(1,11)]
      chr3L = ["chr3."+str(e) for e in range(1,11)]
      chr4L = ["chr4."+str(e) for e in range(1,11)]
      chr5L = ["chr5."+str(e) for e in range(1,11)]
      chr6L = ["chr6."+str(e) for e in range(1,11)]
      chr7L = ["chr7."+str(e) for e in range(1,11)]
      chr8L = ["chr8."+str(e) for e in range(1,11)]
      chr9L = ["chr9."+str(e) for e in range(1,11)]
      chr10L = ["chr10."+str(e) for e in range(1,11)]
      chr11L = ["chr11."+str(e) for e in range(1,11)]
      chr12L = ["chr12."+str(e) for e in range(1,11)]
      chr13L = ["chr13."+str(e) for e in range(1,11)]
      chr14L = ["chr14."+str(e) for e in range(1,11)]
      chr15L = ["chr15."+str(e) for e in range(1,11)]
      chr16L = ["chr16."+str(e) for e in range(1,11)]
      chr17L = ["chr17."+str(e) for e in range(1,11)]
      chr18L = ["chr18."+str(e) for e in range(1,11)]
      chr19L = ["chr19."+str(e) for e in range(1,11)]
      chr20L = ["chr20."+str(e) for e in range(1,11)]
      chr21L = ["chr21."+str(e) for e in range(1,11)]
      chr22L = ["chr22."+str(e) for e in range(1,11)]
      chrXL = ["chrX."+str(e) for e in range(1,11)]

      if(challenge == "train"): chromList = chr2L + chr9L + chr22L
      if(challenge == "test"): chromList = chr1L + chr2L + chr3L + chr4L + chr5L + chr6L + chr7L + chr8L + chr9L + chr10L + chr11L + chr12L + chr13L + chr14L + chr15L + chr16L + chr17L + chr18L + chr19L + chr20L + chr21L + chr22L + chrXL
      if(challenge == "ladder"): chromList = chr1L + chr8L + chr21L

      for chrom in chromList:

        # Parameters
        bedFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/split_"+challenge+"_chrom/"+challenge+"_"+chrom+".bed"
        pfmFileNameList = ",".join(pfmList)
        genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
        footprintBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/footprints/results/"+cell+".bam"
        dnaseBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"+cell+"/"+cell+"_DNase.bam"
        ol = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/"+challenge+"_v5/"
        os.system("mkdir -p "+ol)
        outputFileName = ol+challengeLabel+"."+factor+"."+cell+"."+chrom+".tab"

        # Execution cluster
        myL = "_".join([challenge,factor,cell,chrom])
        clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
        clusterCommand += "-W 120:00 -M 24000 -S 100 -R \"select[hpcwork]\" ./pipeline_v5.zsh "
        clusterCommand += bedFileName+" "+pfmFileNameList+" "+genomeFileName+" "+footprintBamFileName+" "
        clusterCommand += dnaseBamFileName+" "+outputFileName
        os.system(clusterCommand)
        # -P izkf

  tfCellTableFile.close()


