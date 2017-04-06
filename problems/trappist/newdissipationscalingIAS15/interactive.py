import numpy as np
import sys
from subprocess import call
import time

taues = np.logspace(1,3,7)
taues = taues[3:5]
mags = [1.e-3,1.e-2]

for mag in mags:
    for i,taue in enumerate(taues):
        call("python run.py {0} {1} &".format(taue, mag), shell=True)