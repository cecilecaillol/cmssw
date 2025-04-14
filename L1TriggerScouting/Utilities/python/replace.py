import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--skeleton")
parser.add_argument("--input_root")
parser.add_argument("--output_root")
parser.add_argument("--output_py")
input_root = parser.parse_args().input_root
output_root = parser.parse_args().output_root
output_py = parser.parse_args().output_py
skeleton = parser.parse_args().skeleton


fin = open(skeleton, "rt")
#output file to write the result to
fout = open(output_py, "wt")
#for each line in the input file
for line in fin:
   #read replace the string and write to output file
   if '##MYOUTPUT##' in line:
      fout.write(line.replace('##MYOUTPUT##', output_root))
   elif '##MYINPUT##' in line:
      fout.write(line.replace('##MYINPUT##', input_root))
   else:
      fout.write(line)
#close input and output files
fin.close()
fout.close()
