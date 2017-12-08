
# Import
import os
import sys

# Parameters
inLoc = "/hpcwork/izkf/projects/dream_tfbs/local/ChIPseq/peaks/conservative/"
outLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/meme_chip/meme/"
factorList = [
"A549.CEBPB", "A549.CTCF", "A549.E2F6", "A549.GATA3", "A549.MAX", "A549.MYC", "A549.TEAD4", 
"GM12878.ATF2", "GM12878.ATF7", "GM12878.CREB1", "GM12878.E2F1", "GM12878.EGR1", "GM12878.EP300", 
"GM12878.GABPA", "GM12878.MAFK", "GM12878.MAX", "GM12878.RFX5", "GM12878.SPI1", "GM12878.SRF", 
"GM12878.TAF1", "GM12878.TCF12", "GM12878.YY1", "GM12878.ZNF143", "H1-hESC.ATF2", "H1-hESC.ATF3", 
"H1-hESC.CEBPB", "H1-hESC.CREB1", "H1-hESC.CTCF", "H1-hESC.E2F6", "H1-hESC.EGR1", "H1-hESC.EP300", 
"H1-hESC.GABPA", "H1-hESC.MAFK", "H1-hESC.MAX", "H1-hESC.NANOG", "H1-hESC.REST", "H1-hESC.SRF", 
"H1-hESC.TAF1", "H1-hESC.TCF12", "H1-hESC.TEAD4", "H1-hESC.YY1", "H1-hESC.ZNF143", "HCT116.ATF3", 
"HCT116.CEBPB", "HCT116.EGR1", "HCT116.JUND", "HCT116.MAX", "HCT116.SRF", "HCT116.TCF7L2", "HCT116.TEAD4", 
"HCT116.YY1", "HeLa-S3.CEBPB", "HeLa-S3.CTCF", "HeLa-S3.E2F1", "HeLa-S3.E2F6", "HeLa-S3.EP300", 
"HeLa-S3.GABPA", "HeLa-S3.JUND", "HeLa-S3.MAFK", "HeLa-S3.MAX", "HeLa-S3.MYC", "HeLa-S3.REST", 
"HeLa-S3.RFX5", "HeLa-S3.STAT3", "HeLa-S3.TAF1", "HeLa-S3.TCF7L2", "HeLa-S3.ZNF143", "HepG2.ARID3A", 
"HepG2.ATF3", "HepG2.ATF7", "HepG2.CEBPB", "HepG2.CREB1", "HepG2.CTCF", "HepG2.EP300", "HepG2.FOXA1", 
"HepG2.FOXA2", "HepG2.GABPA", "HepG2.HNF4A", "HepG2.JUND", "HepG2.MAFK", "HepG2.MAX", "HepG2.REST", 
"HepG2.SRF", "HepG2.TEAD4", "HepG2.YY1", "HepG2.ZNF143", "IMR-90.CEBPB", "IMR-90.CTCF", "IMR-90.MAFK", 
"K562.ATF3", "K562.ATF7", "K562.CEBPB", "K562.CREB1", "K562.CTCF", "K562.EP300", "K562.JUND", 
"K562.MAX", "K562.MYC", "K562.SRF", "K562.TAF1", "K562.TEAD4", "MCF-7.ATF2", "MCF-7.CTCF", "MCF-7.EGR1", 
"MCF-7.GABPA", "MCF-7.JUND", "MCF-7.MYC", "MCF-7.REST", "MCF-7.RFX5", "MCF-7.TCF12", "Panc1.REST", 
"Panc1.TCF7L2", "SK-N-SH.EP300", "SK-N-SH.GABPA", "SK-N-SH.GATA3", "SK-N-SH.JUND", "SK-N-SH.MAX", 
"SK-N-SH.REST", "SK-N-SH.RFX5", "SK-N-SH.TAF1", "SK-N-SH.TCF12", "SK-N-SH.TEAD4", "SK-N-SH.YY1"
]

# Cell Loop
for factor in factorList:

  # Execution cluster
  #myL = "MEME_"+factor
  #clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
  #clusterCommand += "-W 120:00 -M 48000 -S 100 -P izkf -R \"select[hpcwork]\" ./pipeline.zsh "
  #clusterCommand += inLoc+" "+outLoc+" "+factor
  #os.system(clusterCommand)

  # Execution local
  print "### "+factor+" ------------------------------------------\n"
  sys.stdout.flush()
  clusterCommand = "./pipeline.zsh "+inLoc+" "+outLoc+" "+factor
  os.system(clusterCommand)
  sys.stdout.flush()
  print "----------------------------------------------------------\n\n"
  sys.stdout.flush()


