import sys
import os
import random

# taxon profile creation -- creates taxon profiles to generate simulated reads in MetaSim

def list_files(path):
    # returns a list of names (with extension, without full path) of all files 
    # in folder path
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(str(os.path.join(path, name)))
    return files 

# directory contains reference genomes
indir = sys.argv[1]

# number of dataset to create
ndataset = sys.argv[2]

# name of output taxon profile
outdir = sys.argv[3]

# name of logfile
logfile = "logfile"

files = list_files(indir)

if os.path.isdir(str(outdir)) == False:
    os.mkdir(str(outdir))

#rr is range for random
rr = [1, 5, 10, 20]

for i in range(int(ndataset)):
	output = str(i+1)+'- '+outdir+'.mprf'
	fileOutput = open(output,'w')
	logfile = logfile + str(i+1) +".txt"
	logOut = open(logfile, 'w')
	
	abundance = []
	for x in range(len(files)):
		abundance.append(-1)

	# randomly select 40 genomes to have abundance = 0
	idx = random.sample(range(0, len(files)-1), (len(files)/2))
	for x in range(len(idx)):
		abundance[idx[x]] = 0

	# assign random abundance
	for j in range(len(abundance)):
		if abundance[j] == -1:
				if i <= len(rr)-1:
					rint = random.randint(1, rr[i])
				else:	
					rint = random.randint(1, rr[len(rr)-1])

				abundance[j] = rint
		

	k = 0
	for f in files:
		fi = open(f, "r")

		for line in fi:
		        if line[0] == '>':
		            line = line.strip()
		            part = line.split("|")

		            # create taxon profile by gi
		            fileOutput.write(str(abundance[k]) + " " + part[0].replace(">","") + " " + part[1] + '\n')

		            # create taxon profile by name
		            #name = part[4].split(",")
		            #fileOutput.write(str(random.randint(1,10)) + ' name "' + name[0].strip() + '"\n')

		            # logfile for abundance
		            if abundance[k] != 0:
		            	filename = str(f).split("/")
		            	logOut.write(str(abundance[k]) + " " + str(filename[1]) + '\n')

	    	fi.close()
	    	k = k+1
	fileOutput.close()
	logOut.close()
