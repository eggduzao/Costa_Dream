
# Import
import os
import sys

def maxminFile(fileName):
  tcFile = open(fileName,"r")
  minV = 99999.; maxV = -99999.
  for line in tcFile:
    ll = line.strip().split("\t")
    v = float(ll[3])
    if(v < minV): minV = v
    if(v > maxV): maxV = v
  tcFile.close()
  return minV, maxV

# Input
tcFileName = sys.argv[1]
mpbsFileName = sys.argv[2]
outFileName = sys.argv[3]
pseudocount = 0.00001

# Max and Min
tcMin, tcMax = maxminFile(tcFileName)
mpbsMin, mpbsMax = maxminFile(mpbsFileName)

tcFile = open(tcFileName,"r")
mpbsFile = open(mpbsFileName,"r")
outFile = open(outFileName,"w")
while(True):

  tcl = tcFile.readline()
  mpbsl = mpbsFile.readline()
  if(not tcl): break

  ll = tcl.strip().split("\t")

  tcValue = min(((float(tcl.strip().split("\t")[3])-tcMin)/(tcMax-tcMin))+pseudocount,1) 
  mpbsValue = min(((float(mpbsl.strip().split("\t")[3])-mpbsMin)/(mpbsMax-mpbsMin))+pseudocount,1)

  composedScore = tcValue * mpbsValue
  outFile.write("\t".join([ll[0],ll[1],ll[2],str(composedScore)])+"\n")

tcFile.close()
mpbsFile.close()
outFile.close()


