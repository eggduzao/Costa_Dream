
# Goal: The goal of this script is to retrieve the top X centrimo motifs from meme-chip and convert to our PFM format.

# Import
import os
import sys
from glob import glob

# Input
top = 10
inLocList = glob("/work/eg474423/eg474423_Projects/trunk/Dream/exp/meme_chip/meme/*/")
outLoc = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/motifs/pfm_meme/"

tfList = []
for inLoc in inLocList:

  toRemove = []

  cell = inLoc.split("/")[-2].split(".")[0]
  factor = inLoc.split("/")[-2].split(".")[1]
  if(factor not in tfList): tfList.append(factor)

  # Creating all MEME PFMs
  memeFileName = inLoc+"meme_out/meme.txt"
  memeFile = open(memeFileName,"r")
  resDict = dict()
  count = 0
  for line in memeFile:
    if("position-specific probability matrix" in line):
      nb = line.strip().split(" ")[1]
      outFileName = outLoc+"MEME_"+nb+".pwm"
      toRemove.append(outFileName)
      resDict[outFileName] = [[],[],[],[]]
      count = 3
    elif(count):
      if(count == 3):
        count -= 1
        continue
      elif(count == 2):
        count -= 1
        continue
      elif("--" in line):
        count = 0
        continue
      else:
        ll = [e for e in line.strip().split(" ") if e]
        for i in range(0,4): resDict[outFileName][i].append(ll[i])
  memeFile.close()
  for k in resDict.keys():
    outFile = open(k,"w")
    for v in resDict[k]: outFile.write(" ".join(v)+"\n")
    outFile.close()

  # Creating all DREME PFMs
  dremeFileName = inLoc+"dreme_out/dreme.txt"
  dremeFile = open(dremeFileName,"r")
  resDict = dict()
  flag = False
  for line in dremeFile:
    if("MOTIF " in line):
      nb = line.strip().split(" ")[1]
      outFileName = outLoc+"DREME_"+nb+".pwm"
      toRemove.append(outFileName)
      resDict[outFileName] = [[],[],[],[]]
      count = 3
    elif("letter-probability matrix" in line):
      flag = True
      continue
    elif(flag):
      if(len(line) < 2):
        flag = False
        continue
      ll = [e for e in line.strip().split(" ") if e]
      for i in range(0,4): resDict[outFileName][i].append(ll[i])
  dremeFile.close()
  for k in resDict.keys():
    outFile = open(k,"w")
    for v in resDict[k]: outFile.write(" ".join(v)+"\n")
    outFile.close()

  # Creating all JASPAR PFMs
  jasparFileName = "/work/eg474423/eg474423_Projects/trunk/Dream/exp/meme_chip/JASPAR_CORE_2016_vertebrates.meme"
  jasparFile = open(jasparFileName,"r")
  resDict = dict()
  flag = False
  for line in jasparFile:
    if("MOTIF " in line):
      nb = line.strip().split(" ")[1]
      outFileName = outLoc+"JASPAR_"+nb+".pwm"
      toRemove.append(outFileName)
      resDict[outFileName] = [[],[],[],[]]
      count = 3
    elif("letter-probability matrix" in line):
      flag = True
      continue
    elif(flag):
      if(len(line) < 2):
        flag = False
        continue
      ll = [e for e in "".join(line.strip().split(" ")).split("\t") if e]
      for i in range(0,4): resDict[outFileName][i].append(ll[i])
  jasparFile.close()
  for k in resDict.keys():
    outFile = open(k,"w")
    for v in resDict[k]: outFile.write(" ".join(v)+"\n")
    outFile.close()

  # Iterating on CENTRIMO file and fetching PFMs
  centrimoFileName = inLoc+"centrimo_out/centrimo.txt"
  centrimoFile = open(centrimoFileName,"r")
  centrimoFile.readline(); centrimoFile.readline()
  resTable = []
  for line in centrimoFile:
    ll = [e for e in " ".join(line.strip().split("\t")).split(" ") if e]
    nb = ll[1]
    if(ll[2] != "MEME" and ll[2] != "DREME"): db = "JASPAR"
    else: db = ll[2]
    evalue = int(ll[3].split("e")[1])
    resTable.append([db,nb,evalue])  
  centrimoFile.close()
  resTable = sorted(resTable, key=lambda x: x[2])

  for i in range(0,min(top,len(resTable))):
    outFileN = outLoc+resTable[i][0]+"_"+resTable[i][1]+".pwm"
    outputFileName = outLoc+".".join([cell,factor,str(i+1),resTable[i][0],resTable[i][1],"pwm"])
    os.system("cat "+outFileN+" > "+outputFileName)

  # Removing all pwms
  for e in toRemove: os.system("rm "+e)

# Write output info file name
tfList = sorted(tfList)
outInfoFileName = outLoc+"info.txt"
outInfoFile = open(outInfoFileName,"w")
for tf in tfList:
  allFactorFiles = ",".join([".".join(e.split("/")[-1].split(".")[:-1]) for e in glob(outLoc+"*."+tf+".*.pwm")])
  outInfoFile.write("\t".join([tf,allFactorFiles])+"\n")
outInfoFile.close()


