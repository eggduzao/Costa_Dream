
# Import
import os
import sys

# Input
resLoc = "./results/"
cellList = ["A549", "GM12878", "H1-hESC", "HCT116", "HeLa-S3", "HepG2", "IMR90", 
            "induced_pluripotent_stem_cell", "K562", "liver", "MCF-7", "Panc1", "PC-3", "SK-N-SH"]

# Cell Loop
for cell in cellList:

  print cell
  sys.stdout.flush()
  resFile = resLoc+cell+".bed"
  os.system("awk -F'\\t' 'BEGIN{getline;min=max=$5} NF{ max=(max>$5)?max:$5; min=(min>$5)?$5:min} END{print min,max}' "+resFile)
  sys.stdout.flush()


