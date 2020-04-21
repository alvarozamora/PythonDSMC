import numpy as np


files = 10 #cpus

samples = range(2*10**5, 10**6+2)
samples = range(10**6+2, 2*10**6+2)

samples = range(10)

s = len(samples)

Samples = [ samples[i*s//files:(i+1)*s//files+1] for i in range(files) ]
print(Samples)
args = [[]]
Script =[
"""#!/bin/bash
#SBATCH --job-name=DSMC
#SBATCH --output=DSMCo
#SBATCH --error=DSMCe
#SBATCH --time=24:00:00
#SBATCH -p iric
#SBATCH -N 1
#SBATCH --ntasks-per-node=10

PATH=/home/pizza/yt-conda/bin:$PATH
export PATH
export LD_LIBRARY_PATH=/home/users/kdaegene/hdf5-1.8.20/hdf5/lib:$LD_LIBRARY_PATH
cd /scratch/PI/kipac/pizza/PythonDSMC

ytmod
python DSMC.py """ +str(Samples[f][0]) + ' ' + str(Samples[f][-1]) for f in range(files)]

#print(Script)


for i in range(files):
	filename = "batch"+str(i)


	file = open(filename,'w')

	file.write(Script[i])

	file.close()


