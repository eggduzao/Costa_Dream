
# Import
import os
import sys
import math

def maxminFile(fileName):
  myFile = open(fileName,"r")
  minTc = 99999.; maxTc = -99999.; minBit = 99999.; maxBit = -99999.; minFp = 99999.; maxFp = -99999.
  for line in myFile:
    ll = line.strip().split("\t")
    tc = float(ll[3])
    bit = float(ll[4])
    fp = float(ll[10])
    if(tc < minTc): minTc = tc
    if(tc > maxTc): maxTc = tc
    if(bit < minBit): minBit = bit
    if(bit > maxBit): maxBit = bit
    if(fp < minFp): minFp = fp
    if(fp > maxFp): maxFp = fp
  myFile.close()
  return minTc, maxTc, minBit, maxBit, minFp, maxFp

def sigmoid(x):
  return 1 / (1 + math.exp(-x))

# Input
inFileName = sys.argv[1] 
outFileName = sys.argv[2]

# Parameters
pseudocount = 0.001

# Max and Min
minTc, maxTc, minBit, maxBit, minFp, maxFp = maxminFile(inFileName)

# Iterating in file and calculating score
minScore = 99999.
maxScore = -99999.
tempFileName = outFileName+"_TEMP.bed"
inFile = open(inFileName,"r")
tempFile = open(tempFileName,"w")
for line in inFile:
  ll = line.strip().split("\t")

  #tc = (((float(ll[3]) - minTc) / (maxTc - minTc)) * 20 ) - 10
  #bit = (((float(ll[4]) - minBit) / (maxBit - minBit)) * 20 ) - 10
  tc = (((float(ll[3]) - minTc ) / (maxTc - minTc )) + pseudocount)
  bit = (((float(ll[4]) - minBit ) / (maxBit - minBit )) +pseudocount ) 
  fp = (((float(ll[10]) - minFp ) / (maxFp - minFp ))) + 1 

  score = tc*bit*0.25*fp
  if(score < minScore): minScore = score
  if(score > maxScore): maxScore = score
  tempFile.write("\t".join(ll[:3]+[str(score)])+"\n")
inFile.close()
tempFile.close()

# Correcting score to be within [0,1]
tempFile = open(tempFileName,"r")
outFile = open(outFileName,"w")
for line in tempFile:
  ll = line.strip().split("\t")
  standardizedScore = (float(ll[3]) - minScore) / (maxScore - minScore)
  outFile.write("\t".join(ll[:3]+[str(standardizedScore)])+"\n")
tempFile.close()
outFile.close()
  
# Termination
os.system("rm "+tempFileName)


