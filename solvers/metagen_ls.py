import csv
import sys
import numpy
from numpy import matrix
from sklearn import linear_model
from datetime import datetime
import pylab as pl

# matrix F
infile_F = sys.argv[1]

# vector b
infile_b = sys.argv[2]

outfile = sys.argv[3]

f1 = open(infile_F, "r")
reader = csv.reader(f1)

# g is number of genomes
g = 0
li = f1.readline()
g = len(li.split(","))

# n is number of k-mers
n = 0
for row in reader:
    n += 1
n = n + 1

f1.close()

f = open(infile_F, "r")
s1 = datetime.now()

print "Reading matrix F: %s" %infile_F
F = numpy.zeros(shape=(n,g))
name = []

i = 0
for line in f.readlines():
    line = line.strip()
    part = line.split(",")
    if (part[0] != "K-mer"):
    	for j in range(1,len(part)):
    		F[i][j-1] = float(part[j])
    	i = i + 1
    else:
        for idx in range(1,len(part)):
            name.append(part[idx])
    

s2 = datetime.now()

f2 = open(infile_b, "r")
print "Reading vector b: %s" %infile_b
b = numpy.zeros(shape=(n,1))

i = 0
for line in f2.readlines():
	line = line.strip()
	part = line.split(",")
	if (part[0] != "K-mer"):
		b[i][0] = (float(part[len(part)-1]))
		i = i + 1

clf = linear_model.LinearRegression()
clf.__init__(fit_intercept=False, normalize=False, copy_X=True)
clf.fit(F, b)
r = clf.coef_

out = open(outfile,'w')
if (len(r[0]) > 0):
	for i in range(0, g-1):
		out.write("%s, %f\n" %(name[i], r[0,i]))
	out.close()
