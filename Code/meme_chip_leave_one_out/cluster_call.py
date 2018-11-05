
# Import
import os
import sys
from math import floor

# Input
chipHalfExt = "50"
tfCellFileName = "/work/eg474423/Costa_Dream/exp/tf_cell_table.txt"
memeDBFileName = "/work/eg474423/Costa_Dream/exp/meme_chip/JASPAR_CORE_2016_vertebrates.meme"
genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
tfPeakLoc = "/hpcwork/izkf/projects/dream_tfbs/local/ChIPseq/peaks/conservative/"
outLoc = "/work/eg474423/Costa_Dream/exp/meme_chip_leave_one_out/results/"

# Fetching cell list for each TF
tfCellDict = dict()
tfCellFile = open(tfCellFileName, "r")
for line in tfCellFile:
  ll = line.strip().split("\t")
  tfCellDict[ll[0]] = ll[1].split(",")
tfCellFile.close()
factorList = sorted(tfCellDict.keys())

# Iterating on each factor
for factor in factorList:

  for cell in tfCellDict[factor]:

    inCellList = ",".join([x for x in tfCellDict[factor] if x != cell])

    # Execution cluster
    myL = "MEMELOO_"+cell+"_"+factor
    clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
    clusterCommand += "-W 120:00 -M 48000 -S 100 -R \"select[hpcwork]\" ./pipeline.zsh "
    clusterCommand += factor+" "+cell+" "+chipHalfExt+" "+memeDBFileName+" "+genomeFileName+" "+tfPeakLoc+" "+outLoc+" "+inCellList
    os.system(clusterCommand)
    # -P izkf

    #print "### "+factor+" / "+cell+" ------------------------------------------\n"
    #sys.stdout.flush()
    # Command
    #sys.stdout.flush()
    #print "----------------------------------------------------------\n\n"
    #sys.stdout.flush()


