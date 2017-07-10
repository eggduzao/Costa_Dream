
# Import
import os
import sys

# Input

mpbsFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/test_regions.blacklistfiltered.bed"
bamLoc = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"
outputLoc = "/work/ig440396/project/Dream/exp/tc/results/"
cellList = ["A549", "GM12878", "H1-hESC", "HCT116", "HeLa-S3", "HepG2", "IMR90", 
            "induced_pluripotent_stem_cell", "K562", "liver", "MCF-7", "Panc1", "PC-3", "SK-N-SH"]

# Cell Loop
for cell in cellList:

  bamFileName = bamLoc+cell+"/"+cell+"_DNase.bam"
  outputFileName = outputLoc+cell+".bed"

  # Execution
  myL = "TC_"+cell
  clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
  clusterCommand += "-W 120:00 -M 12000 -S 100 -P izkf -R \"select[hpcwork]\" ./pipeline.zsh "
  clusterCommand += mpbsFileName+" "+bamFileName+" "+outputFileName
  os.system(clusterCommand)


