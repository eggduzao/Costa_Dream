
# Import
import os
import sys

# Input
inLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"
outLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/kmer_bias/bias/"
cellList = ["A549", "GM12878", "H1-hESC", "HCT116", "HeLa-S3", "HepG2", "IMR90", 
            "induced_pluripotent_stem_cell", "K562", "liver", "MCF-7", "Panc1", "PC-3", "SK-N-SH"]

# Cell Loop
for cell in cellList:

  # Execution
  myL = "KMER_"+cell
  clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
  clusterCommand += "-W 10:00 -M 24000 -S 100 -P izkf -R \"select[hpcwork]\" ./pipeline.zsh "
  clusterCommand += inLoc+" "+outLoc+" "+cell
  os.system(clusterCommand)


