#from __future__ import division
#from math import *
#import numpy as np
#from ofblockmeshdicthelper import *
#import fileinput
#import sys
#import os
#import matplotlib.pyplot as plt


# ========================================================================================
f=open('./cases/airfoil_run/postProcessing/forceCoeffs/0/forceCoeffs.dat','r')
lines=f.readlines()
split=lines[-1].split()
Cd=float(split[2])
Cl=float(split[3])
f.close()
CdCl=abs(Cd/Cl)
# ========================================================================================
f = open('results.txt', 'w')
f.write(str(CdCl))
#f.write("\n")
f.close()

