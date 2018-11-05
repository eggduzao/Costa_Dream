
import os
import sys
from glob import glob
from Bio import motifs

inList = glob("./*.EP300.*.pwm")
for pwmFileName in inList:
  pwmFile = open(pwmFileName,"r")
  logo_file_name = pwmFileName[:-3]+"png"
  pwm = motifs.read(pwmFile, "pfm")
  pwm.weblogo(logo_file_name, format="png_print", stack_width = "medium", color_scheme = "color_classic")
  pwmFile.close()


