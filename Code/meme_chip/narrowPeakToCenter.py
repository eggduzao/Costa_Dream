import os
import sys
halfExt = int(sys.argv[1])
inFile = open(sys.argv[2],"r")
outFile = open(sys.argv[3],"w")
for line in inFile:
  ll = line.strip().split("\t")
  mid = int(ll[1]) + int(ll[9])
  outFile.write("\t".join([ll[0],str(mid-halfExt),str(mid+halfExt)])+"\n")
inFile.close()
outFile.close()
