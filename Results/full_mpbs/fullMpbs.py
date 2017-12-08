
# Import
import os
import sys
from Bio import motifs
from pysam import Fastafile
from MOODS import search

# Motif class
class Motif:
  def __init__(self, input_file_name):

    # Standardize input file to be only the nucleotide frequencies

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

# Input
bedFileName = sys.argv[1]
pfmFileNameList = sys.argv[2].split(",")
genomeFileName = sys.argv[3]
outputFileName = sys.argv[4]

# Fetching motifs
motifList = []
for pfmFileName in pfmFileNameList:
  motifList.append(Motif(pfmFileName))
globalMin = min([e.min for e in motifList])

bedFile = open(bedFileName,"r")
outFile = open(outputFileName,"w")
genomeFile = Fastafile(genomeFileName)
for line in bedFile:

  ll = line.strip().split("\t")
  chrName = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])

  try: sequence = str(genomeFile.fetch(chrName, p1, p2))
  except Exception:
    print "Exception 1 raised in "+line
    outFile.write("\t".join([ll[0],ll[1],ll[2],str(globalMin)])+"\n")
    continue

  maxValue = globalMin

  try:
    for res in search(sequence, [e.pssm_list for e in motifList], [e.min for e in motifList], both_strands=True):
      for (position, score) in res:
        if(score > maxValue): maxValue = score
  except Exception:
    print "Exception 2 raised in "+line
    outFile.write("\t".join([ll[0],ll[1],ll[2],str(globalMin)])+"\n")
    continue
  
  outFile.write("\t".join([ll[0],ll[1],ll[2],str(maxValue)])+"\n")

bedFile.close()
outFile.close()
genomeFile.close()


