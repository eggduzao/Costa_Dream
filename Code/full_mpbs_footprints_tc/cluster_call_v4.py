
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

    ########### Should be REMOVED later TODO XXX #####
    #if((challenge == "train") and (factor != "CTCF" and factor != "FOXA1" and factor != "FOXA2" and factor != "HNF4A" and factor != "TAF1")): continue
    #if((challenge == "test") and (factor != "CTCF" and factor != "EGR1" and factor != "FOXA1" and factor != "FOXA2" and factor != "GABPA" and factor != "HNF4A" and factor != "JUND" and factor != "MAX" and factor != "NANOG" and factor != "REST" and factor != "TAF1")): continue
    #if((challenge == "ladder") and (factor != "CTCF" and factor != "EGR1")): continue
    ##################################################

    trainCellList = ll[1].split(",")
    testCellList = ll[challengeCol].split(",")

    # Creating PFM list based on training cells
    pfmLoc = "/work/eg474423/Costa_Dream/exp/motifs/pfm_hocomoco_positive_negative/"
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
      if(challenge == "train"): chromList = ["chr2", "chr9", "chr22"]
      if(challenge == "test"): chromList = [e.split("/")[-1].split("_")[-1].split(".")[0] for e in glob(bedLoc+"*_chr*.bed")]
      if(challenge == "ladder"): chromList = ["chr1", "chr8", "chr21"]

      ########### Should be REMOVED later TODO XXX #####
      #chr1L = ["chr1."+str(e) for e in range(1,11)]
      #chr2L = ["chr2."+str(e) for e in range(1,11)]
      #chr3L = ["chr3."+str(e) for e in range(1,11)]
      #chr4L = ["chr4."+str(e) for e in range(1,11)]
      #chr5L = ["chr5."+str(e) for e in range(1,11)]
      #chr6L = ["chr6."+str(e) for e in range(1,11)]
      #chr7L = ["chr7."+str(e) for e in range(1,11)]
      #chr8L = ["chr8."+str(e) for e in range(1,11)]
      #chr10L = ["chr10."+str(e) for e in range(1,11)]
      #chr11L = ["chr11."+str(e) for e in range(1,11)]
      #chr12L = ["chr12."+str(e) for e in range(1,11)]
      #chr13L = ["chr13."+str(e) for e in range(1,11)]
      #chrXL = ["chrX."+str(e) for e in range(1,11)]

      #if(challenge == "train"): chromList = chr2L
      #if(challenge == "test"):
      #  if(factor == "CTCF"): chromList = chr1L + chr2L + chr3L + chr4L
      #  if(factor == "EGR1"): chromList = chr1L + chr2L
      #  if(factor == "FOXA1"): chromList = chr1L + chr2L + chr3L + chr4L + chr5L + chr6L + chr7L + chr8L + chrXL
      #  if(factor == "FOXA2"): chromList = chr1L + chr2L + chr3L + chr4L + chr5L + chr6L + chr7L + chrXL
      #  if(factor == "GABPA"): chromList = chr1L + chr2L
      #  if(factor == "HNF4A"): chromList = chr1L + chr2L + chr3L + chr4L + chr5L + chr6L + chr7L + chrXL + chr10L + chr11L + chr12L + chr13L
      #  if(factor == "JUND"): chromList = chr1L + chr2L + chr3L + chr4L + chr5L
      #  if(factor == "MAX"): chromList = chr1L + chr2L + chr3L + chr6L + chr7L
      #  if(factor == "NANOG"): chromList = chr1L + chr2L
      #  if(factor == "REST"): chromList = chr1L + chr2L + chr3L
      #  if(factor == "TAF1"): chromList = chr1L + chr2L + chr3L + chr4L + chr5L + chr6L + chr8L + chrXL
      #if(challenge == "ladder"): chromList = chr1L
      ##################################################

      for chrom in chromList:

        # Parameters
        bedFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/split_"+challenge+"_chrom/"+challenge+"_"+chrom+".bed"
        pfmFileNameList = ",".join(pfmList)
        genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
        footprintBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/footprints/results/"+cell+".bam"
        dnaseBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"+cell+"/"+cell+"_DNase.bam"
        ol = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/"+challenge+"_v4.2/"
        os.system("mkdir -p "+ol)
        outputFileName = ol+challengeLabel+"."+factor+"."+cell+"."+chrom+".tab"

        # Execution cluster
        myL = "_".join([challenge,factor,cell,chrom])
        clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
        clusterCommand += "-W 120:00 -M 24000 -S 100 -R \"select[hpcwork]\" ./pipeline_4.zsh "
        clusterCommand += bedFileName+" "+pfmFileNameList+" "+genomeFileName+" "+footprintBamFileName+" "
        clusterCommand += dnaseBamFileName+" "+outputFileName
        os.system(clusterCommand)
        # -P izkf

  tfCellTableFile.close()


