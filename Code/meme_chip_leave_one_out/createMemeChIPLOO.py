
# Import
import os
import sys
from math import floor

# Input
factor = sys.argv[1]
cell = sys.argv[2]
chipHalfExt = int(sys.argv[3])
memeDBFileName = sys.argv[4]
genomeFileName = sys.argv[5]
tfPeakLoc = sys.argv[6]
outLoc = sys.argv[7]
inCellList = sys.argv[8].split(",")

# Initialization
if(len(inCellList) <= 0): sys.exit(0)
totalPeakEach = int(floor(500.0/len(inCellList)))
totalPeaks = totalPeakEach * len(inCellList)
newOL = outLoc+cell+"."+factor+"/"
os.system("mkdir -p "+newOL)
toRemove = []

# Get TF peaks
chipFileNameA = newOL+"chipA.bed"; toRemove.append(chipFileNameA)
chipFileNameB = newOL+"chipB.bed"; toRemove.append(chipFileNameB)
chipFileNameC = newOL+"chipC.bed"; toRemove.append(chipFileNameC)
for c in inCellList:
  chipFileName = tfPeakLoc+"ChIPseq."+c+"."+factor+".conservative.train.narrowPeak.gz"
  os.system("gunzip -c "+chipFileName+" > "+chipFileNameA)
  inFile = open(chipFileNameA,"r")
  outFile = open(chipFileNameB,"w")
  for line in inFile:
    ll = line.strip().split("\t")
    mid = int(ll[1]) + int(ll[9])
    outFile.write("\t".join([ll[0],str(mid-chipHalfExt),str(mid+chipHalfExt)])+"\n")
  inFile.close()
  outFile.close()
  os.system("shuf -n "+str(totalPeakEach)+" "+chipFileNameB+" >> "+chipFileNameC)

# Get fasta file
faFileName = newOL+"myfa.bed"; toRemove.append(faFileName)
os.system("fastaFromBed -fi "+genomeFileName+" -bed "+chipFileNameC+" -fo "+faFileName)

# Apply MEME
os.system("meme-chip -oc "+newOL+" -db "+memeDBFileName+" -nmeme "+str(totalPeaks)+" -meme-mod zoops -meme-minw 5 -meme-maxw 15 -meme-nmotifs 5 "+faFileName)
# -meme-p 4

# Termination
for e in toRemove: os.system("rm "+e)


