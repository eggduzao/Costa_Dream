
import os
import sys

myL = "test_dream"
clusterCommand = "bsub -J "+myL+" -o "+myL+"_out.txt -e "+myL+"_err.txt "
clusterCommand += "-W 200:00 -M 24000 -S 100 -R \"select[hpcwork]\" -P izkf ./pipecall.zsh"
os.system(clusterCommand)

