
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
tcHalfWindowVec = [50]
maxReadPile = 1000

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

# Function calculate tag count
def tag_count(chrName, start, end, bamFile, tcHalfWindow=None):
  total = 0
  if(tcHalfWindow):
    mid = (start+end)/2
    p1exttc = max(mid - tcHalfWindow,0)
    p2exttc = mid + tcHalfWindow
  else:
    p1exttc = start
    p2exttc = end
  for read in bamFile.fetch(reference=chrName, start=p1exttc, end=p2exttc):
    total += 1
    if(total >= maxReadPile): break
  return total

# Function to calculate overlap
def overlap(t1, t2):
  return max(0, min(t1[1], t2[1]) - max(t1[0], t2[0]))

# Function to write output
def writeOutput(ll,regionTagCountVec,resVec,outFile):
  finalResVec = []
  for vec in resVec:
    for e in vec: finalResVec.append(e)
  outFile.write("\t".join([ll[0],ll[1],ll[2]]+[str(e) for e in regionTagCountVec]+[str(e) for e in finalResVec])+"\n")

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

# Making header naming
metricForEachMotif = ["BITSCORE", "BITSCORE_FOOTPRINT", "PROTECTION", "RELATIVE_MOTIF_POSITION", "OVERLAP_FOOTPRINT"]
pfmFileNameListShort = [e.split("/")[-1].split(".")[0] for e in pfmFileNameList]
header = "TC"
counterPos = 1
counterNeg = 1
for e in pfmFileNameListShort:
  myvec = []
  if("POS." in e):
    factorname = "\tDENOVO_POS_"+str(counterPos)
    counterPos += 1
  elif("NEG." in e):
    factorname = "\tDENOVO_NEG_"+str(counterNeg)
    counterNeg += 1
  else: factorname = "\tHOCOMOCO"
  for k in metricForEachMotif: myvec.append(factorname+"_"+k)
  header += "".join(myvec)

# Iterating on main bed file
bedFile = open(bedFileName,"r")
outFile = open(outputFileName,"w")
outFile.write(header+"\n")
genomeFile = Fastafile(genomeFileName)

for line in bedFile:

  # Fetching line
  ll = line.strip().split("\t")
  chrName = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])

  # Starting result structures
  regionTagCount = 0
  resVec = []
  for m in motifList:
    vec = []
    # For each motif print:
    #1. Bitscore.
    #2. Bitscore Footprint.
    #3. Protection.
    #4. Relative motif position.
    #5. Amount of overlap between FP and motif.
    vec.append(globalMin); vec.append("NA"); vec.append("NA"); vec.append("NA"); vec.append("NA")
    resVec.append(vec)
  
  # Evaluating Overall TC
  regionTagCountVec = ["NA"]
  try:
    for i in range(0,len(tcHalfWindowVec)):
      tcHalfWindow = tcHalfWindowVec[i]
      regionTagCountVec[i] = tag_count(chrName, p1, p2, dnaseBam, tcHalfWindow)
  except Exception: 
    print "Exception TC raised in "+line
    writeOutput(ll,regionTagCountVec,resVec,outFile)
    continue

  # Fetching sequence
  try: sequence = str(genomeFile.fetch(chrName, p1, p2))
  except Exception:
    print "Exception SEQUENCE raised in "+line
    writeOutput(ll,regionTagCountVec,resVec,outFile)
    continue

  # Performing motif matching
  for i in range(0,len(motifList)):

    m = motifList[i]
    for res in search(sequence, [m.pssm_list], [m.min], absolute_threshold=True, both_strands=True):
      for (position, score) in res:
        if(score > resVec[i][0]):
          resVec[i][0] = score
          resVec[i][3] = abs(position)

    # Fetching absolute position of best motif match
    mp1 = p1 + position
    mp2 = p1 + m.len

    # Evaluating PS
    p1_ext = mp1 - m.len; p2_ext = mp2 + m.len
    if(p1_ext < 0): continue
    nl = tag_count(chrName, p1_ext, mp1, dnaseBam)
    nc = tag_count(chrName, mp1, mp2, dnaseBam)
    nr = tag_count(chrName, mp2, p2_ext, dnaseBam)
    resVec[i][2] = round(((nr+nl)/2.0)-nc, 6)

  # Fetching footprints
  try:
    footprints = fpBam.fetch(reference=chrName, start=p1, end=p2)
    for fp in footprints:
      sequence = str(genomeFile.fetch(chrName, fp.pos-5, fp.aend+5))

      for i in range(0,len(motifList)):
        m = motifList[i]
        flag = False
        for res in search(sequence, [m.pssm_list], [m.min], absolute_threshold=True, both_strands=True):
          for (position, score) in res:
            if(score > resVec[i][1]):
              resVec[i][1] = score
              flag = True

        if(flag):

          # Fetching absolute position of best motif match
          mp1 = p1 + position
          mp2 = p1 + m.len

          # Evaluating overlap
          resVec[i][4] = overlap([mp1, mp2], [fp.pos, fp.aend])

  except Exception:
    print "Exception FOOTPRINTS raised in "+line
    writeOutput(ll,regionTagCountVec,resVec,outFile)
    continue

  # Writing results
  writeOutput(ll,regionTagCountVec,resVec,outFile)

# Termination
bedFile.close()
outFile.close()
genomeFile.close()
dnaseBam.close()
fpBam.close()


