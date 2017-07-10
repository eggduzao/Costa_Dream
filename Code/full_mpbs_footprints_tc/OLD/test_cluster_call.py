
# Import
import os
import sys
from glob import glob

# Iterating on TF-CELL table
tfCellTableFileName = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/tf_cell_table.txt"
tfCellTableFile = open(tfCellTableFileName,"r")
tfCellTableFile.readline()
for line in tfCellTableFile:

  # Line parsing
  ll = line.strip().split("\t")
  if(ll[3] == "."): continue
  factor = ll[0]
  trainCellList = ll[1].split(",")
  testCellList = ll[3].split(",")

  # Creating PFM list based on training cells
  pfmLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/motifs/pfm_meme/"
  pfmInfoFileName = pfmLoc+"info.txt"
  pfmInfoFile = open(pfmInfoFileName,"r")
  pfmList = []
  for line2 in pfmInfoFile:
    ll2 = line2.strip().split("\t")

  pfmInfoFile.close()

  for tcell in trainCellList:
    for e in glob(pfmLoc+tcell+"."+factor+".*.pwm"): pfmList.append(e)

  # Iterating on test cells
  for cell in testCellList:

    # Iterating on chromosomes
    chromList = ["chr"+str(e) for e in range(1,23)+["X"]]
    for chrom in chromList:

      # Parameters
      bedFileName = "/hpcwork/izkf/projects/dream_tfbs/local/annotations/split_test_chrom/test_"+chrom+".bed"
      pfmFileNameList = ",".join(pfmList)
      genomeFileName = "/hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.fa"
      footprintBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/footprints/results/"+cell+".bam"
      dnaseBamFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"+cell+"/"+cell+"_DNase.bam"
      outputFileName = "/hpcwork/izkf/projects/dream_tfbs/exp/full_mpbs_footprints_tc/test/F."+factor+"."+cell+"."+chrom+".tab"

      if(chrom in ["chr"+str(e) for e in range(1,8)+["X"]]):
        proj = "-P izkf"
        time = "300:00"
      else:
        proj = ""
        time = "120:00"

      # Execution cluster
      myL = "FINAL_"+factor+"_"+cell+"_"+chrom
      clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
      clusterCommand += "-W "+time+" -M 12000 -S 100 "+proj+" -R \"select[hpcwork]\" ./pipeline.zsh "
      clusterCommand += bedFileName+" "+pfmFileNameList+" "+genomeFileName+" "+footprintBamFileName+" "
      clusterCommand += dnaseBamFileName+" "+outputFileName
      os.system(clusterCommand)

tfCellTableFile.close()


