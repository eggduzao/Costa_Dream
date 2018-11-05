
import os
import sys

cell = "K562"
factor = "ARID3A"
loc = "/home/egg/eg474423_Projects/trunk/Dream/exp/motifs/pfm_meme/"
pfmList = [loc+"HepG2.ARID3A.2.MEME.1.pwm", loc+"HepG2.ARID3A.3.DREME.AGRKGGC.pwm"]

#bedFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/mytest2.bed"
bedFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"+cell+"/DNase_Peaks.bed"
pfmFileNameList = ",".join(pfmList)
genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
footprintBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/footprints/results/"+cell+".bam"
dnaseBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"+cell+"/"+cell+"_DNase.bam"
outputFileName = "./test.tab"

# Execution
os.system("python fullMpbs.py "+bedFileName+" "+pfmFileNameList+" "+genomeFileName+" "+footprintBamFileName+" "+dnaseBamFileName+" "+outputFileName)


