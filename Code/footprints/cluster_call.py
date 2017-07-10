
# Import
import os
import sys

# Input
matrixLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/footprints/matrix/"
biasLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/kmer_bias/bias/table_format/"
outLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/footprints/results/"
cellList = ["A549", "GM12878", "H1-hESC", "HCT116", "HeLa-S3", "HepG2", "IMR90", 
            "induced_pluripotent_stem_cell", "K562", "liver", "MCF-7", "Panc1", "PC-3", "SK-N-SH"]

# Cell Loop
for cell in cellList:

  #peakFile = "/hpcwork/izkf/projects/dream_tfbs/exp/dnase/"+cell+"/DNase_Peaks.bed"
  #print cell
  #sys.stdout.flush()
  #for i in range(1,23)+["X"]:
  #  os.system("grep chr"+str(i)+" "+peakFile+" | head -1 ")

  # Execution
  #myL = "HINT_"+cell
  #clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
  #clusterCommand += "-W 6:00 -M 24000 -S 100 -P izkf -R \"select[hpcwork]\" ./pipeline.zsh "
  #clusterCommand += matrixLoc+" "+biasLoc+" "+outLoc+" "+cell
  #os.system(clusterCommand)

  # Execution local
  print "### "+cell+" ------------------------------------------\n"
  sys.stdout.flush()
  clusterCommand = "./pipeline.zsh "+matrixLoc+" "+biasLoc+" "+outLoc+" "+cell
  os.system(clusterCommand)
  sys.stdout.flush()
  print "----------------------------------------------------------\n\n"
  sys.stdout.flush()


