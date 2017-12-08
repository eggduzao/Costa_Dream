
# Import
import os
import sys
from glob import glob

# Iterating on TF-CELL table
tfCellTableFileName = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/tf_cell_table.txt"
tfCellTableFile = open(tfCellTableFileName,"r")
tfCellTableFile.readline()
for line in tfCellTableFile:

  # Line parsing
  ll = line.strip().split("\t")
  if(ll[2] == "."): continue
  factor = ll[0]
  trainCellList = ll[1].split(",")
  testCellList = ll[2].split(",")

  # Creating PFM list based on training cells
  pfmLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/motifs/pfm_meme/"
  pfmList = []
  for tcell in trainCellList:
    for e in glob(pfmLoc+tcell+"."+factor+".*.pwm"): pfmList.append(e)

  # Iterating on test cells
  for cell in testCellList:

    # Parameters
    bedFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/ladder_regions.blacklistfiltered.bed"
    pfmFileNameList = ",".join(pfmList)
    genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
    footprintBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/footprints/results/"+cell+".bam"
    dnaseBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"+cell+"/"+cell+"_DNase.bam"
    outputFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/ladder/L."+factor+"."+cell+".tab"

    # Execution cluster
    myL = "BS_"+factor+"_"+cell
    clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
    clusterCommand += "-W 300:00 -M 12000 -S 100 -P izkf -R \"select[hpcwork]\" ./pipeline.zsh "
    clusterCommand += bedFileName+" "+pfmFileNameList+" "+genomeFileName+" "+footprintBamFileName+" "
    clusterCommand += dnaseBamFileName+" "+outputFileName
    os.system(clusterCommand)

tfCellTableFile.close()


