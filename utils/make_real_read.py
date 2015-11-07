import csv
import sys
import os

infile = sys.argv[1]
outfile = sys.argv[2]

f = open(infile, "r")

fileOutput = open(outfile,'w')

for line in f:
    if line[0] == '>':
            line = line.strip()
            fileOutput.write(line.replace('>','@') + '\n')
    else:
            line = line.strip()
            len_line = len(line)
            content = '!' * len_line
            fileOutput.write(line + '\n')
            fileOutput.write('+\n')
            fileOutput.write(content + '\n')
            
f.close()
fileOutput.close()
