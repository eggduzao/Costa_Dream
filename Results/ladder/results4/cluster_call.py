
# Import
import os
import sys
from glob import glob

# Input
inLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/ladder_foot/"
outLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/ladder/results4/"

# Loop
inFileList = glob(inLoc+"*.tab")
#inFileList = [ inLoc+"L.ARID3A.K562.tab" ]
for inF in inFileList:
  print inF
  # TF and CELL
  tf = inF.split("/")[-1].split(".")[1]
  cell = inF.split("/")[-1].split(".")[2]

  # Parameters
  inFileName = inF
  outFileName = outLoc+"L."+tf+"."+cell+".tab"

  # Execution cluster
  myL = "LR3_"+tf+"_"+cell
  clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
  clusterCommand += "-W 120:00 -M 12000 -S 100 -R \"select[hpcwork]\" -P izkf ./pipeline.zsh "
  clusterCommand += inFileName+" "+outFileName
  os.system(clusterCommand)
  print(clusterCommand)
  # -P izkf

# /hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/ladder/L.ARID3A.K562.tab
# /hpcwork/izkf/projects/dream_tfbs/exp/ladder/results2/L.ARID3A.K562.tab
#python ladder3.py /hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/ladder/L.ARID3A.K562.tab /hpcwork/izkf/projects/dream_tfbs/exp/ladder/results3/L.ARID3A.K562.tab
#gzip L.ARID3A.K562.tab
#scp eg474423@$CLUSTER:/hpcwork/izkf/projects/dream_tfbs/exp/ladder/results3/L.ARID3A.K562.tab.gz /home/egg/Desktop/

