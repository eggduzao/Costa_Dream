
# Import
import os
import sys
import math
from copy import deepcopy
import numpy
from bedtools import Interval, IntervalFile
from Bio import motifs
from pysam import Samfile, Fastafile
from MOODS import search

# Input
originalDNaseLocation = "/hpcwork/izkf/projects/dream_tfbs/local/essential_training_data/DNASE/bam/" # Change this line
genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa" # Change this line
annotationBedFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/split_test_chrom/test_chr22.bed" # Change this line
dreamFlag = "F"
cell = "HepG2"
factor = "ATF2"

# Fix paths
if(originalDNaseLocation[-1] != "/"): originalDNaseLocation+="/"

# Parameters
processedDNaseLocation = os.path.abspath("./dnase/")
#footprintLocation = os.path.abspath("./footprints/")
footprintLocation = os.path.abspath("/hpcwork/izkf/projects/dream_tfbs/exp/footprints/results/")
motifLocation = os.path.abspath("./motifs/")
outLoc = os.path.abspath("./out/")

if(processedDNaseLocation[-1] != "/"): processedDNaseLocation+="/"
if(footprintLocation[-1] != "/"): footprintLocation+="/"
if(motifLocation[-1] != "/"): motifLocation+="/"
if(outLoc[-1] != "/"): outLoc+="/"

qValueThresh = 10.0
dnaseBamFileName = processedDNaseLocation+cell+"/"+cell+"_DNase.bam"
hsFileName = processedDNaseLocation+cell+"/DNase_Peaks.bed"

def readPfmList(factor, inFileName):
  inFile = open(inFileName,"r")
  res = []
  for line in inFile:
    ll = line.strip().split("\t")
    if(ll[0] == factor):
      res = deepcopy(ll)
      break
  inFile.close()
  return res

motifListFileNameFpr = motifLocation+"info_fpr.txt"
motifListFileName = motifLocation+"info.txt"
motifListVecFpr = readPfmList(factor, motifListFileNameFpr)
motifListVec = readPfmList(factor, motifListFileName)
fprCutoff = motifListVecFpr[2]
pfmFileNameListFpr = [motifLocation+e+".pwm" for e in motifListVecFpr[1].split(",")]
pfmFileNameList = [motifLocation+e+".pwm" for e in motifListVec[1].split(",")]

tcHalfWindow = 50
maxReadPile = 200

toRemove = []

##########################################################################
# DNase pre-process
##########################################################################

# Pre-process dnase and call peaks
os.system("./dnase_preprocess.sh "+originalDNaseLocation+" "+processedDNaseLocation+" "+cell)

# Merging overlapping peaks and keeping max score
inFileName = processedDNaseLocation+cell+"/Peaks/Peaks_peaks.narrowPeak"
mergFileName = processedDNaseLocation+cell+"/Peaks/merged.bed"
toRemove.append(mergFileName)
os.system("mergeBed -c 9 -o max -i "+inFileName+" > "+mergFileName)

# Removing q-value < thresh
inFile = open(mergFileName,"r")
outFile = open(hsFileName,"w")
counter = 1
for line in inFile:
  ll = line.strip().split("\t")
  if(ll[0] == "chrY" or ll[0] == "chrM" or "random" in ll[0] or "Un" in ll[0] or "_" in ll[0]): continue
  if(float(ll[3]) > qValueThresh):
    outFile.write("\t".join([ll[0],ll[1],ll[2],"p"+str(counter),ll[3]])+"\n")
    counter += 1
inFile.close()
outFile.close()

##########################################################################
# Evaluate Footprints
##########################################################################

# We provide the evaluated footprints in "footprintLocation" since its calculation is time-consuming. Please contact the authors for all the scripts for the generation of bias-correction matrices and HINT-BC application.

footprintFileName = footprintLocation+cell+".bed"
footprintBamFileName = footprintLocation+cell+".bam"

##########################################################################
# Evaluate bit-score inside footprints
##########################################################################

bitInsideFootFileName = outLoc+"bitInsideFoot.bed"
toRemove.append(bitInsideFootFileName)

# Motif class
class Motif:
  def __init__(self, input_file_name):
    input_file = open(input_file_name,"r")
    self.pfm = motifs.read(input_file, "pfm")
    self.pwm = self.pfm.counts.normalize(0.1)
    input_file.close()
    self.len = len(self.pfm)
    background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
    self.pssm = self.pwm.log_odds(background)
    self.pssm_list = [self.pssm[e] for e in ["A","C","G","T"]]
    self.max = self.pssm.max
    self.min = self.pssm.min

# Fetching motifs
motifList = []
for pfmFileName in pfmFileNameListFpr:
  motifList.append(Motif(pfmFileName))
globalMin = min([e.min for e in motifList])
motifSize = motifList[0].len
bedFile = open(footprintFileName,"r")
outFile = open(bitInsideFootFileName,"w")
genomeFile = Fastafile(genomeFileName)
for line in bedFile:
  ll = line.strip().split("\t")
  chrName = ll[0]; p1 = int(ll[1]) + motifSize; p2 = int(ll[2]) + motifSize
  try: sequence = str(genomeFile.fetch(chrName, p1, p2))
  except Exception:
    print "Exception 1 raised in "+line
    outFile.write("\t".join([ll[0],ll[1],ll[2],str(globalMin)])+"\n")
    continue
  for res in search(sequence, [e.pssm_list for e in motifList], fprCutoff, threshold_from_p=False,both_strands=True,convert_log_odds=False,log_base=2):
      if len(res)==0:
        continue
      maxScore=-1000
      for (position, score) in res:
        maxScore=max(maxScore,score)
      outFile.write("\t".join([ll[0],ll[1],ll[2],str(maxScore)])+"\n")
bedFile.close()
outFile.close()
genomeFile.close()

##########################################################################
# Extract Features (TC, bit-score and footprint)
##########################################################################

featuresFileName = outLoc+"featuresFile.bed"
toRemove.append(featuresFileName)

#############
# Functions / Classes
#############

# Motif class
class Motif:
  def __init__(self, input_file_name):
    input_file = open(input_file_name,"r")
    self.pfm = motifs.read(input_file, "pfm")
    self.pwm = self.pfm.counts.normalize(0.0001)
    input_file.close()
    self.len = len(self.pfm)
    background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
    self.pssm = self.pwm.log_odds(background)
    self.pssm_list = [self.pssm[e] for e in ["A","C","G","T"]]
    self.max = self.pssm.max
    self.min = self.pssm.min

# Function to read pileup
def tag_count(chrName, start, end, bamFile, tcHalfWindow):
  total = 0
  mid = (start+end)/2
  p1exttc = max(mid - tcHalfWindow,0)
  p2exttc = mid + tcHalfWindow
  for read in bamFile.fetch(reference=chrName, start=p1exttc, end=p2exttc):
    total += 1
    if(total >= maxReadPile): break
  return total

# Function to calculate overlap
def overlap(t1, t2):
  return max(0, min(t1[1], t2[1]) - max(t1[0], t2[0]))

# Function to write output
def writeOutput(ll,regionTagCount,resVec,outFile):
  outFile.write("\t".join([ll[0],ll[1],ll[2],str(regionTagCount)]+[str(e) for e in resVec])+"\n")

#############
# Execution
#############

# Fetching motifs & evaluating global minimum
motifList = []
for pfmFileName in pfmFileNameList:
  motifList.append(Motif(pfmFileName))
globalMin = min([e.min for e in motifList])

# Opening BAM files for DNase-seq and footprints (footprint bed files were converted to BAM for efficiency)
dnaseBam = Samfile(dnaseBamFileName,"rb")
fpBam = Samfile(footprintBamFileName,"rb")

# Iterating on main bed file
bedFile = open(annotationBedFileName,"r")
outFile = open(featuresFileName,"w")
genomeFile = Fastafile(genomeFileName)
for line in bedFile:

  # Fetching line
  ll = line.strip().split("\t")
  chrName = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])

  # Starting result structures
  regionTagCount = 0
  resVec = [globalMin,0,0,0,0,0] # BIT-SCORE, MOTIF_P1, MOTIF_P2, FP_OVERLAP, FP_P1, FP_P2
  counter = 0

  # Evaluating Overall TC
  try: regionTagCount = tag_count(chrName, p1, p2, dnaseBam, tcHalfWindow)
  except Exception: 
    print "Exception TC raised in "+line
    writeOutput(ll,regionTagCount,resVec,outFile)
    continue

  # Fetching sequence
  try: sequence = str(genomeFile.fetch(chrName, p1, p2))
  except Exception:
    print "Exception SEQUENCE raised in "+line
    writeOutput(ll,regionTagCount,resVec,outFile)
    continue

  # Fetching footprints
  try: footprints = fpBam.fetch(reference=chrName, start=p1, end=p2)
  except Exception:
    print "Exception FOOTPRINTS raised in "+line
    writeOutput(ll,regionTagCount,resVec,outFile)
    continue

  # Best mpbs
  maxPos = -99999
  maxValue = globalMin
  maxMotifLen = -1

  # Performing motif matching and footprint overlapping
  for res in search(sequence, [e.pssm_list for e in motifList], [e.min for e in motifList], absolute_threshold=True, both_strands=True):

    for (position, score) in res:

      if(score > maxValue):
        maxValue = score
        maxPos = position
        maxMotifLen = motifList[counter].len

    counter += 1
    if(counter == len(resVec)): break

  # Overlap and Writting output
  if(maxValue != globalMin):
    # t1 = MPBS
    t1 = [0,0]
    t2Write = [0,0]
    if(maxPos >= 0): t1 = [p1+maxPos, p1+maxPos+maxMotifLen]
    else: t1 = [p1-maxPos, p1-maxPos+maxMotifLen]
    maxOverlap = 0
    for f in footprints:
      # t2 = footprint
      t2 = [f.pos,f.aend]
      overlapN = overlap(t1, t2)
      if(overlapN > maxOverlap):
        maxOverlap = overlapN
        t2Write[0] = t2[0]
        t2Write[1] = t2[1]
    resVec = [maxValue, t1[0], t1[1], maxOverlap, t2Write[0], t2Write[1]]
    writeOutput(ll,regionTagCount,resVec,outFile)
  else: writeOutput(ll,regionTagCount,resVec,outFile)

# Termination
bedFile.close()
outFile.close()
genomeFile.close()
dnaseBam.close()
fpBam.close()

##########################################################################
# Add footprint column to features file
##########################################################################

newFeaturesFileName = outLoc+"newFeaturesFile.bed"
toRemove.append(newFeaturesFileName)

min_value=1000000
inFile = open(bitInsideFootFileName,"r")
for line in inFile:
  line=line.strip("\n")
  ll = line.split("\t")
  min_value=min(min_value,float(ll[3])-0.01)
inFile.close()

footprints = IntervalFile(bitInsideFootFileName)
inFile = open(featuresFileName,"r")
outFile = open(newFeaturesFileName,"w")
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

##########################################################################
# Integrating all metrics (TC, bit-score and footprint)
##########################################################################

outFileName = outLoc+dreamFlag+"."+factor+"."+cell+".tab"

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

# Parameters
pseudocount = 0.001

# Max and Min
minTc, maxTc, minBit, maxBit, minFp, maxFp = maxminFile(newFeaturesFileName)

# Iterating in file and calculating score
minScore = 99999.
maxScore = -99999.
tempFileName = outFileName+"_TEMP.bed"
toRemove.append(tempFileName)
inFile = open(newFeaturesFileName,"r")
tempFile = open(tempFileName,"w")
for line in inFile:
  ll = line.strip().split("\t")
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

##########################################################################
# Termination
##########################################################################

os.system("gzip -c "+outFileName+" > "+outFileName+".gz")
#for e in toRemove: os.system("rm "+e)


