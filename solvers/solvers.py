# LP, L1, L2 solvers without reducing
import csv
import sys
from gurobipy import *
from datetime import datetime
import scipy.optimize
import numpy
from numpy import matrix
from numpy import array,linspace,zeros,eye,concatenate,sum as SUM,linalg
from scipy.optimize import nnls
from sklearn import linear_model
import pylab as pl

# matrix F
infile1 = sys.argv[1]
# b vector
infile2 = sys.argv[2]
# output file for LP
outfile_lp = sys.argv[3]
# output file for L1
outfile_l1 = sys.argv[4]
# output file for L2
outfile_ls = sys.argv[5]

print "Read F matrix: %s" %infile1
f1 = open(infile1, "r")

F = {}
name = []

# n is number of k-mers
s1 = datetime.now()
n = 0
for line1 in f1.readlines():
    line1 = line1.strip()
    #print(line)
    part1 = line1.split(",")
    #print part
    if (part1[0] != "K-mer"):
        for j in range(1,len(part1)):
        	F[n,j-1] = float(part1[j])
        n = n + 1
    else:
        for idx in range(1,len(part1)):
            name.append(part1[idx])

# g is number of genomes
g = len(part1)-1

print("Number of K-mer = %d\n" %n)
print("Number of Genomes = %d\n" %g)
#for row in range(n):
#    for col in range(g):
#        print(F[row, col],)
        #print("  ")
#    print("\n")

s2 = datetime.now()

print "Read b vector: %s" %infile2
f2 = open(infile2, "r")

b = []

for line2 in f2.readlines():
    line2 = line2.strip()
    #print(line)
    part2 = line2.split(",")
    #print part
    if (part2[0] != "K-mer"):
        b.append(float(part2[len(part2)-1])) 

#for i in range(n):
#    print b[i]


# solve by LP
try:
    print "\nBuilding LP model......"
    s3 = datetime.now()
    m1 = Model("metagen_lp")

    # Variables
    x1 = [None] * (g)
    y1 = [None] * (n)
    for i in range(g):
        x1[i] = m1.addVar(lb = 0, vtype = GRB.CONTINUOUS, name=name[i])
    for j in range(n):
        y1[j] = m1.addVar(lb = 0, vtype = GRB.CONTINUOUS, name="y_lp"+str(j))
    
    m1.update()

    # Objective
    obj = LinExpr() 
    for i in range(n):
        obj += y1[i]
    
    m1.setObjective(obj, GRB.MINIMIZE)

    # Constraints
    for i in range(n):
        constr = LinExpr()
        for j in range(g):
            constr += F[i,j] * x1[j]
        m1.addConstr(constr + y1[i] == b[i])
        
    print "Finish building model"
    s4 = datetime.now()
    print s4 - s3

    s5 = datetime.now()
    print "\n"
    m1.optimize()
    m1.write("file1.lp")
        
    s6 = datetime.now()

    if m1.SolCount > 0:
        out_lp = open(outfile_lp,'w')
        info = ""
        for v in m1.getVars():    
            out_lp.write('%s, %f \n' %(v.varName, float(v.x)))
        out_lp.close()

except GurobiError as e:
    print('Error reported')
    print e.errno

"""
# solve by L1
try:
    print "\nBuilding LP model for L1......"
    s3 = datetime.now()
    m2 = Model("lpfroml1")

    # Variables
    x2 = [None] * (g)
    y2 = [None] * (n)
    for i in range(g):
        x2[i] = m2.addVar(lb = 0, vtype = GRB.CONTINUOUS, name=name[i])
    for j in range(n):
        y2[j] = m2.addVar(lb = 0, vtype = GRB.CONTINUOUS, name="y_l1"+str(j))
    
    m2.update()

    # Objective
    obj = LinExpr() 
    for i in range(n):
        obj += y2[i]
    
    m2.setObjective(obj, GRB.MINIMIZE)

    # Constraints
    nc = 0
    for i in range(n):
            constr1 = LinExpr()
            constr2 = LinExpr()
            for j in range(g):
                constr1 += F[i,j] * x2[j]
                constr2 += -1.0 * F[i,j] * x2[j]
            nc += 2
            m2.addConstr(constr1 - y2[i] <= b[i])
            m2.addConstr(constr2 - y2[i] <= -1.0 * b[i])
        
    print "Finish building model"
    s4 = datetime.now()
    print s4 - s3

    print "Number of constr = %d\n" %nc

    s5 = datetime.now()
    print "\n"
    m2.optimize()
    m2.write("file2.lp")
              

    s6 = datetime.now()

    if m2.SolCount > 0:
        out_l1 = open(outfile_l1,'w')
        info = ""
        for v in m2.getVars():    
            out_l1.write('%s, %f \n' %(v.varName, float(v.x)))
        out_l1.close()

except GurobiError as e:
    print('Error reported')
    print e.errno
"""

# solve by L2
print "\nRunning Non-negative Least Square solver...."
F3 = numpy.zeros(shape=(n,g))
#b3 = numpy.zeros(shape=(n,1))

for i in range(0, n):
    for j in range(0, g):
        F3[i][j] = F[i, j]

sol = scipy.optimize.nnls(F3, b)
#weights = normalise(sol[0])

out = open(outfile_ls,'w')
for i in range(0, len(sol[0])):
    out.write("%s, %.30f\n" %(name[i], round(sol[0][i], 30)))
out.close()
