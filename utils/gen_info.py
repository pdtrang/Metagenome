import csv
import sys
import os

# get genome information (gi: genome index) from the first line

def list_files(path):
    # returns a list of names (with extension, without full path) of all files 
    # in folder path
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(str(os.path.join(path, name)))
    return files 

indir = sys.argv[1]
outfile = sys.argv[2]

files = list_files(indir)

fileOutput = open(outfile,'w')

for f in files:
		fi = open(f, "r")
		filename = str(f).split("/")

		for line in fi:
		        if line[0] == '>':
		            line = line.strip()
		            part = line.split("|")

		            # create taxon profile by gi
		            fileOutput.write(str(filename[len(filename)-1]) + ", " + str(part[1]) + '\n')

	    	fi.close()

fileOutput.close()
