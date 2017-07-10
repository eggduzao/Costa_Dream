
# Import
import os
import sys
from glob import glob

# Input
tfcelltableFileName = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/tf_cell_table.txt"
tcLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/ladder/tc/results/"
mpbsLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs/hocomoco/"
outLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/ladder/results1/"

# Loop
tfcelltableFile = open(tfcelltableFileName,"r")
tfcelltableFile.readline()
for line in tfcelltableFile:

  ll = line.strip().split("\t")
  if(ll[2] == "."): continue
  tf = ll[0]
  for cell in ll[2].split(","):
  
    # Parameters
    tcFileName = tcLoc+cell+".bed"
    mpbsFileName = mpbsLoc+tf+".bed"
    outFileName = outLoc+"L."+tf+"."+cell+".tab"

    # File exists
    if(not os.path.exists(tcFileName)):
      print "("+tf+","+cell+") = "+tcFileName+" does not exist"+"\n"
      continue
    if(not os.path.exists(mpbsFileName)):
      print "("+tf+","+cell+") = "+mpbsFileName+" does not exist"+"\n"
      continue

    # Execution cluster
    myL = "LR1_"+tf+"_"+cell
    clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
    clusterCommand += "-W 120:00 -M 12000 -S 100 -P izkf -R \"select[hpcwork]\" ./pipeline.zsh "
    clusterCommand += tcFileName+" "+mpbsFileName+" "+outFileName
    os.system(clusterCommand)

tfcelltableFile.close()


