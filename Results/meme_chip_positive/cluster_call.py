
# Import
import os
import sys
from math import floor

# Input
chipHalfExt = "50"
totalPeaks = "500"
tfCellFileName = "/work/eg474423/Costa_Dream/exp/tf_cell_table.txt"
memeDBFileName = "/work/eg474423/Costa_Dream/exp/meme_chip/JASPAR_CORE_2016_vertebrates.meme"
genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
tfPeakLoc = "/hpcwork/izkf/projects/dream_tfbs/local/ChIPseq/peaks/conservative/"
outLoc = "/work/eg474423/Costa_Dream/exp/meme_chip_positive/results/"

# Fetching cell list for each TF
tfCellDict = dict()
tfCellFile = open(tfCellFileName, "r")
tfCellFile.readline()
for line in tfCellFile:
  ll = line.strip().split("\t")
  #if(ll[0] != "EGR1" and ll[0] != "CTCF"): continue # Apply only for EGR1 and CTCF
  #if(ll[0] == "EGR1" or ll[0] == "CTCF"): continue # Do not apply for EGR1 and CTCF
  tfCellDict[ll[0]] = ll[1].split(",")
tfCellFile.close()
factorList = sorted(tfCellDict.keys())

# Iterating on each factor
for factor in factorList:

  cellList = ",".join(tfCellDict[factor])

  # This code below is to call the script in the cluster
  myL = "MEMPOS_"+factor
  clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
  clusterCommand += "-W 120:00 -M 48000 -S 100 -R \"select[hpcwork]\" ./pipeline.zsh "
  clusterCommand += factor+" "+cellList+" "+chipHalfExt+" "+totalPeaks+" "+memeDBFileName+" "+genomeFileName+" "+tfPeakLoc+" "+outLoc
  os.system(clusterCommand)
  ## Option not used because I submit to the regular cluster: -P izkf

  # This code below is to call the script locally
  #print "### "+factor+" / "+cellList+" ------------------------------------------\n"
  #sys.stdout.flush()
  #clusterCommand = "./pipeline.zsh "+factor+" "+cellList+" "+chipHalfExt+" "+totalPeaks+" "+memeDBFileName+" "+genomeFileName+" "+tfPeakLoc+" "+outLoc
  #os.system(clusterCommand)
  #sys.stdout.flush()
  #print "----------------------------------------------------------\n\n"
  #sys.stdout.flush()


