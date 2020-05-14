from numpy import loadtxt
import numpy as np
from sys import argv
a=loadtxt(argv[1])
MAE=np.average(np.abs(a[:,2:]),axis=0)
error = np.average(MAE)
print("#",error)
for value in MAE:
    print(value)
