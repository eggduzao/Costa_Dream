
# Import
import os
import sys

# Input
factor = sys.argv[1]
cellList = sys.argv[2].split(",")
chipHalfExt = int(sys.argv[3])
dnaseHalfExt = int(sys.argv[4])
totalPeaks = int(sys.argv[5])
memeDBFileName = sys.argv[6]
genomeFileName = sys.argv[7]
tfPeakLoc = sys.argv[8]
dnasePeakLoc = sys.argv[9]
outLoc = sys.argv[10]

# Initialization
newOL = outLoc+factor+"/"
os.system("mkdir -p "+newOL)
toRemove = []

chipFileNameA = newOL+"chipA.bed"; toRemove.append(chipFileNameA)
chipFileNameB = newOL+"chipB.bed"; toRemove.append(chipFileNameB)
chipFileNameC = newOL+"chipC.bed"; toRemove.append(chipFileNameC)

dnaseFileNameA = newOL+"dnaseA.bed"; toRemove.append(dnaseFileNameA)

intFileNameA = newOL+"intA.bed"; toRemove.append(intFileNameA)
intFileNameB = newOL+"intB.bed"; toRemove.append(intFileNameB)
intFileNameC = newOL+"intC.bed"; toRemove.append(intFileNameC)
outFileTT = open(intFileNameC,"w")
outFileTT.close()

faFileName = newOL+"myfa.fa"; toRemove.append(faFileName)

# Get TF peaks
for cell in cellList:

  chipFileName = tfPeakLoc+"ChIPseq."+cell+"."+factor+".conservative.train.narrowPeak.gz"
  os.system("gunzip -c "+chipFileName+" > "+chipFileNameA)
  inFile = open(chipFileNameA,"r")
  outFile = open(chipFileNameB,"w")
  for line in inFile:
    ll = line.strip().split("\t")
    mid = int(ll[1]) + int(ll[9])
    outFile.write("\t".join([ll[0],str(mid-chipHalfExt),str(mid+chipHalfExt)])+"\n")
  inFile.close()
  outFile.close()
  os.system("sort -k1,1 -k2,2n "+chipFileNameB+" > "+chipFileNameC)

  # Get DNase peaks
  dnaseFileName = dnasePeakLoc+cell+"/DNase_Peaks.bed"

  os.system("sort -k1,1 -k2,2n "+dnaseFileName+" > "+dnaseFileNameA)

  # Perform subtraction

  os.system("intersectBed -a "+dnaseFileNameA+" -b "+chipFileNameC+" -wa -v > "+intFileNameA)
  inFile = open(intFileNameA,"r")
  outFile = open(intFileNameB,"w")
  for line in inFile:
    ll = line.strip().split("\t")
    mid = (int(ll[1]) + int(ll[2])) / 2
    outFile.write("\t".join([ll[0],str(mid-dnaseHalfExt),str(mid+dnaseHalfExt)])+"\n")
  inFile.close()
  outFile.close()
  os.system("shuf -n "+str(int(totalPeaks/len(cellList))+1)+" "+intFileNameB+" >> "+intFileNameC)

  # Get fasta file
  os.system("fastaFromBed -fi "+genomeFileName+" -bed "+intFileNameC+" -fo "+faFileName)

# Apply MEME
os.system("meme-chip -oc "+newOL+" -db "+memeDBFileName+" -nmeme "+str(totalPeaks)+" -meme-mod zoops -meme-minw 5 -meme-maxw 15 -meme-nmotifs 5 "+faFileName)
# -meme-p 4

# Termination
for e in toRemove: os.system("rm "+e)


