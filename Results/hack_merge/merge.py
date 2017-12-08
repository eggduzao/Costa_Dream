
# Import
import os
import sys
import math
from bedtools import Interval, IntervalFile

# Input
inFileName = sys.argv[1] 
inFileName2 = sys.argv[2] 
outFileName = sys.argv[3]

min_value=1000000
inFile = open(inFileName2,"r")
for line in inFile:
  line=line.strip("\n")
  ll = line.split("\t")
  min_value=min(min_value,float(ll[3])-0.01)
inFile.close()

footprints = IntervalFile(inFileName2)
inFile = open(inFileName,"r")
outFile = open(outFileName,"w")
for line in inFile:
  line=line.strip("\n")
  ll = line.split("\t")

  chr=ll[0]
  pos1=int(ll[1])
  pos2=int(ll[2])

  query = Interval(chr,pos1,pos2)
  
  score = min_value
  for h in footprints.search(query):
    score=max(score,float(h.name))
  
  outFile.write(line+"\t"+str(score)+"\n")

outFile.close()
  


