import glob
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sod
import pdb
import yt
yt.enable_parallelism()


#Sod Exact, Riemann
pl = 1.01
pr = 1.0
Pl = 1
Pr = 0.99 #8/10.*pl/pr
gamma = 5./3
dustFrac = 0.0
npts = 100
t = 0.1
print(t)
left_state = (Pl, pl,0)
right_state = (Pr, pr, 0.)

positions, regions, values = sod.solve(left_state=left_state, 
    right_state=right_state, geometry=(0., 1., 0.5), t=t, 
    gamma=gamma, npts=npts, dustFrac=dustFrac)


#DSMC
#files = glob.glob("SodData/rho*")
#print(len(files))
#pdb.set_trace()
boxes = 100
hist = np.zeros(boxes)

actual = 0
n = 2*10**5
while n < 2*10**6+2:
	n = str(n)
	while len(n) < 8:
		n = '0' + n
	try:
		h1 = np.load("SodData/rho"+str(n)+".npy")
		hist += h1
		actual += 1
	except:
		actual += 0
		print(int(n)-2*10**5, actual)
	finally:
		n = int(n) + 1

hist /= actual

print("Actual Files = ", actual)

#Plotting
plt.figure()
plt.plot(values['x'], hist,'b', label="DSMC")
plt.plot(values['x'], values['rho'], 'k', label='Riemann')
plt.xlabel('x')
plt.ylabel('Density')
plt.grid()
plt.savefig("Density.png")

#plt.figure()
#plt.plot(values['x'], hist,'b', label="DSMC")
#plt.plot(values['x'], values['rho'], 'k', label='Riemann')
#plt.xlabel('x')
#plt.ylabel('Density')
#plt.grid()
#plt.savefig("Density.png")
