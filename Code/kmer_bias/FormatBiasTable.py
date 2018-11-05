
import os
import sys
from glob import glob

inFileList = glob("/work/eg474423/eg474423_Projects/trunk/Dream/exp/kmer_bias/bias/table/*.txt")
outLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/kmer_bias/bias/table_format/"

for inFileName in inFileList:

  if("_All" in inFileName): continue

  inName = inFileName.split("/")[-1]
  outFileName = outLoc + inName

  inFile = open(inFileName,"r")
  outFile = open(outFileName,"w")
  inFile.readline()
  for line in inFile:
    ll = line.strip().split("\t")
    outFile.write("\t".join(ll[:2])+"\n")
  inFile.close()
  outFile.close()


