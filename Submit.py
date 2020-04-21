import os

files = 10
for i in range(files):
	os.system('sbatch Batches3/batch'+str(i))
