
# Import
import os
import sys
from glob import glob

# Input
inLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/test/"
footprintMotif = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_tc/hocomoco/"
outLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/test_foot/"

# Loop
inFileList = glob(inLoc+"*.tab")
for inF in inFileList:
  # TF and CELL
  tf = inF.split("/")[-1].split(".")[1]
  cell = inF.split("/")[-1].split(".")[2]

  # Parameters
  inFileName = inF
  inFileName2 = footprintMotif+tf+"."+cell+".bed"
  outFileName = outLoc+"L."+tf+"."+cell+".tab"

  # Execution cluster
  myL = "LR3_"+tf+"_"+cell
  clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
  clusterCommand += "-W 120:00 -M 12000 -S 100 -R \"select[hpcwork]\" -P izkf ./pipeline.zsh "
  clusterCommand += inFileName+" "+inFileName2+" "+outFileName
  #os.system(clusterCommand)
  print(clusterCommand)
  # -P izkf


