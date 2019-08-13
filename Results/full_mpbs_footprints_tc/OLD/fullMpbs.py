
###########################################################
# Import
###########################################################

# Import
import os
import sys
from Bio import motifs
from pysam import Samfile
from pysam import Fastafile
from MOODS import search

###########################################################
# Input
###########################################################

# Input
bedFileName = sys.argv[1]
pfmFileNameList = sys.argv[2].split(",")
genomeFileName = sys.argv[3]
footprintBamFileName = sys.argv[4]
dnaseBamFileName = sys.argv[5]
outputFileName = sys.argv[6]

# Parameters
tcHalfWindow = 50

###########################################################
# Functions / Classes
###########################################################

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
  for read in bamFile.fetch(reference=chrName, start=p1exttc, end=p2exttc): total += 1
  return total

# Function to calculate overlap
def overlap(t1, t2):
  return max(0, min(t1[1], t2[1]) - max(t1[0], t2[0]))

# Function to write output
def writeOutput(ll,regionTagCount,resVec,outFile):
  outFile.write("\t".join([ll[0],ll[1],ll[2],str(regionTagCount)]))
  for vec in resVec: outFile.write("\t"+"\t".join([str(e) for e in vec]))
  outFile.write("\n")

###########################################################
# Execution
###########################################################

# Fetching motifs & evaluating global minimum
motifList = []
for pfmFileName in pfmFileNameList:
  motifList.append(Motif(pfmFileName))
globalMin = min([e.min for e in motifList])

# Opening BAM files for DNase-seq and footprints (footprint bed files were converted to BAM for efficiency)
dnaseBam = Samfile(dnaseBamFileName,"rb")
fpBam = Samfile(footprintBamFileName,"rb")

# Iterating on main bed file
bedFile = open(bedFileName,"r")
outFile = open(outputFileName,"w")
genomeFile = Fastafile(genomeFileName)
for line in bedFile:

  # Fetching line
  ll = line.strip().split("\t")
  chrName = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])

  # Starting result structures
  regionTagCount = 0
  resVec = [[globalMin,0,0] for e in motifList]
  counter = 0

  # Evaluating Overall TC
  try: regionTagCount = tag_count(chrName, p1, p2, dnaseBam, tcHalfWindow) #
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

  # Performing motif matching and footprint overlapping
  for res in search(sequence, [e.pssm_list for e in motifList], [e.min for e in motifList], absolute_threshold=True, both_strands=True):

    maxValueOverlap = globalMin
    maxValueNoOverlap = globalMin

    for (position, score) in res:

      for f in footprints:
        # t1 = MPBS
        if(position >= 0): t1 = [p1+position,p1+position+motifList[counter].len]
        else: t1 = [p1-position,p1-position+motifList[counter].len]
        # t2 = footprint
        t2 = [f.pos,f.aend]
        overlapN = overlap(t1, t2)
        if((overlapN > 0) and (score > maxValueOverlap)):
          maxValueOverlap = score
          maxValueOverlapTC = tag_count(chrName, position, position+motifList[counter].len, dnaseBam, tcHalfWindow)
            
      if(score > maxValueNoOverlap): maxValueNoOverlap = score

    if((maxValueOverlap > 0) and (overlapN >= motifList[counter].len/2)):
      resVec[counter][0] = maxValueOverlap # MPBS max score
      resVec[counter][1] = overlapN # MPBS+FP overlap
      resVec[counter][2] = maxValueOverlapTC # overlapping FP TC
    else:
      resVec[counter][0] = maxValueNoOverlap # MPBS max score
      resVec[counter][1] = 0 # MPBS+FP overlap
      resVec[counter][2] = 0 # overlapping FP TC

    counter += 1
    if(counter == len(resVec)): break

  # Writing output
  writeOutput(ll,regionTagCount,resVec,outFile)

# Termination
bedFile.close()
outFile.close()
genomeFile.close()
dnaseBam.close()
fpBam.close()


