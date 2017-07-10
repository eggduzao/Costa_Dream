
import os
import sys
from glob import glob

toSplit = 10.0

def file_len(fname):
  i = 0
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i + 1

inFileList = glob("./*.bed")

for inFileName in inFileList:

  flen = file_len(inFileName)
  each = int(flen / toSplit)

  inFile = open(inFileName,"r")
  outFile = None
  counter = 0
  counter2 = 1
  for line in inFile:

    if((counter%each == 0) and (counter2 <= (int(toSplit)))): 
      if(outFile): outFile.close()
      outFile = open(inFileName[:-3]+str(counter2)+".bed","w")
      counter2 += 1
      counter = 0

    outFile.write(line)
    counter += 1
  
  inFile.close()
  if(outFile): outFile.close()


