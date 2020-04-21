import numpy as np
#import matplotlib.pyplot as plt
import pdb
import copy
import time
import multiprocessing as mp
#import sodexact.sod as sod
from functools import partial
mp.cpu_count()


# In[2]:


def Advect(P, deltat, dX):
    #print("Advecting")
    dtsim = np.zeros(len(P))
    DT = np.ones(len(P))*deltat
    
    iteration = 0
    while not np.all(np.isclose(dtsim,DT)):
        #print( (dtsim==DT).sum())
        #Ensure Particles are in box
      
        assert((P[:,:3] >= 0).sum() + (P[:,:3] <= dX).sum() == 6*len(P) )
                  
        dT = DT-dtsim
        dT[dT<0] = 0. #Machine Error
        
        iteration += 1
        #print("Iteration = ", iteration)#, "dT = ", dT)
        #Left Side Collision
        dt1 = -P[:,:3]/P[:,3:6]
        dt1[dt1 <= 0] = 2*deltat #Don't want negative values (past) or machine error 
        R1 = np.argmin(dt1,axis=1)
        dT1 = np.min(dt1,axis=1)
        
        
        #Right Side Collision
        dt2 = (dX[None,:]-P[:,:3])/P[:,3:6]
        dt2[dt2 <= 0] = 2*deltat #Don't want negative values (past) or machine error 
        R2 = np.argmin(dt2,axis=1)
        dT2 = np.min(dt2,axis=1)
        
        R = np.argmin(np.array([dT, dT1, dT2]), axis=0)
        dT = np.min(np.array([dT, dT1, dT2]), axis=0)
       
        #Advect
        #print("Advecting")
        #print(P[:,:3],dT[:,None]*P[:,3:],P[:,:3]+dT[:,None]*P[:,3:])
        
        #Advect
        P[:, :3] = P[:, :3] + dT[:,None]*P[:, 3:6]
        P[R == 1, R1[R == 1]] = 0.0
        P[R == 2, R2[R == 2]] = dX[R2[R==2]]
        
        #P[:,:3] = P[:,:3] + dT[:,None]*P[:,3:6]
        
        #Reflect if necessary
        P[R == 1, 3+R1[R == 1]] = -P[R == 1, 3+R1[R == 1]]
        P[R == 2, 3+R2[R == 2]] = -P[R == 2, 3+R2[R == 2]]
        
        dtsim += dT
    
    return P


# In[3]:


def Sort(P,boxes):
    #print("Sorting")
    bins = np.arange(boxes)/boxes
    inds = np.digitize(P[:,2], bins)-1
    
    Boxes = [[] for box in range(boxes)]
    for i in range(len(P)):
        Boxes[inds[i]].append(P[i])
        
    return Boxes
        
def Vmax(Box):
    vmax = np.zeros(3)
    vmin = np.zeros(3)
    for i in range(len(Box)):
        if i == 0:
            vmax = Box[i][3:6]
            vmin = Box[i][3:6]
        else:
            vmax = np.max([vmax, Box[i][3:6]],axis=0)
            vmin = np.min([vmin, Box[i][3:6]],axis=0)
            
    v2 = (vmax-vmin)**2
    v = np.sqrt(v2.sum())
    return v


def ParallelCollide(Boxes, dX, dT, s, Ne):
    #print("Colliding")
    p = mp.Pool()
    func = partial(Collide,dX,dT,s,Ne)
    boxes = p.map(func, [[Box] for Box in Boxes])
    p.close()
    p.join()
    return [box[0] for box in boxes]


def Collide(dX,dT,s,Ne,Boxes):
    bx = 0
    for Box in range(len(Boxes)):
        bx += 1
        Nc = len(Boxes[Box])
        Vc = np.product(dX)/boxes
        
        
        nc = Nc*Ne/Vc
        vmax = Vmax(Boxes[Box])
        #sigma = s/MFP(s,n)**2/boxes/boxes
        sigma = s
        dV = 0
        Ncand = Nc*Nc*sigma*dT*vmax*Ne/(2*(Vc-dV))*Y
        #if bx == 1:
        #print("Nc = ", Nc, "Candidates = ",Ncand)
        ncand = 0
        
        while ncand < Ncand and Nc > 1:
            ncand += 1
        
            i = int(np.random.uniform()*Nc)
            j = int(np.random.uniform()*Nc)
            while j == i:
                j = int(np.random.uniform()*Nc)
            
            vrel = np.sqrt(((Boxes[Box][i][3:6]-Boxes[Box][j][3:6])**2).sum())
            
            coin = np.random.uniform()
            if vrel/vmax > coin:
                a = 2*np.pi*np.random.uniform()
                b = 2*np.random.uniform()-1
                
                v1 = copy.deepcopy(Boxes[Box][i][3:6])
                v2 = copy.deepcopy(Boxes[Box][j][3:6])
                
                Boxes[Box][i][3:6] = (v1 + v2)/2 + vrel*np.array([np.cos(a)*np.sqrt(1-b**2),np.sin(a)*np.sqrt(1-b**2),b])/2
                Boxes[Box][j][3:6] = (v1 + v2)/2 - vrel*np.array([np.cos(a)*np.sqrt(1-b**2),np.sin(a)*np.sqrt(1-b**2),b])/2
                
                if np.isclose(v1,Boxes[Box][i][3:6]).sum()==3:
                    print("No Momentum Exchange")
                
    return Boxes


def MFP(n,s):
    return 1/(np.sqrt(2)*s*n)

def MFT(n,s,T):
    return 1/np.sqrt(8/np.pi*T)/boxes

def HMS(T):
    H = int(T//3600)
    M = int(T//60)
    S = int(T - H*3600 - M*60)
    
    return str(H)+" hours, "+str(M)+" minutes, "+str(S)+" seconds"
     
def MakeParticles(N, p, rho, z1, z2):
    P = np.ones((N,6))
    P[:,0] = np.random.uniform(size=N)*dX[0]
    P[:,1] = np.random.uniform(size=N)*dX[1]
    P[:,2] = (np.random.uniform(size=N)*(z2-z1) + z1)*dX[2]

    P[:,3] = np.random.normal(size=N)*np.sqrt(p/rho)#*KnEff/Kn
    P[:,4] = np.random.normal(size=N)*np.sqrt(p/rho)#*KnEff/Kn
    P[:,5] = np.random.normal(size=N)*np.sqrt(p/rho)#*KnEff/Kn
    #
    return [np.array(p) for p in P]


# In[4]:


#Shock ICs 0.01
pl = 1.0
pr = 0.1
Pl = 1.0
Pr = 0.1 #8/10.*pl/pr
g = 5/3
M = pl*0.5 + pr*0.5

#Maximal Knudsen of Base Case
boxes = 10000
bins = 100
lam = 1.0/boxes
Kn = 2*np.abs(pl-pr)/(pl+pr)*bins/boxes

#Initialization
n = 1e27
s = 1/(np.sqrt(2)*n*lam)
b2 = 2/3*np.pi*np.sqrt(s/np.pi)**3
Y = (1 + 0.05556782*b2*n + 0.01394451*b2**2*n**2 - 0.0013396*b2**3*n**3)/(1 - 0.56943218*b2*n + 0.08289011*b2**2*n**2)



print("Knudsen = ",Kn, ", s = ", s, ", Y = ",Y)

#box size (not grid size), made into cubes
dx = dy = 1
dz = boxes
dX = np.array([dx,dy,dz])/boxes


dt = MFT(n,s,Pl/pl)/Y
#print("KnEff = ",KnEff)

#Ep = 0.0005
#Ep = 0.002
Ep =  0.01

amp = 20
simcell = 100*amp
samples = int(1.0/Ep**2/simcell)+1

def Init(Simcell):
    #Particle Initialization
    
    left  = [ MakeParticles(10*amp, Pl, pl, box/boxes, (box+1)/boxes) for box in range(boxes//2)]
    right = [ MakeParticles(1*amp, Pr, pr, box/boxes+0.5, (box+1)/boxes+0.5) for box in range(boxes//2)]


    Boxes = left + right
    P = np.concatenate(Boxes,axis=0)
    N = len(P)
    Ne = n*MFP(s,n)**3*boxes**3*np.product(dX)/N
    #print("N = ", N, " , Ne = ", Ne, ", boxes = ", boxes, ", Samples = ", np.round(samples,1))


    '''
    init = np.random.uniform(size=N)
    init = (init < pl/(pl+pr)).astype(int).sum()

    P = np.ones((N,6))
    P[:,0] = np.random.uniform(size=N)*dx
    P[:,1] = np.random.uniform(size=N)*dy
    P[:init,2] = np.random.uniform(size=init)*0.5
    P[init:,2] = (np.random.uniform(size=N-init)*0.5 + 0.5)

    P[:init,3] = np.random.normal(size=init)*np.sqrt(Pl/pl)*KnEff/Kn
    P[init:,3] = np.random.normal(size=N-init)*np.sqrt(Pr/pr)*KnEff/Kn
    P[:init,4] = np.random.normal(size=init)*np.sqrt(Pl/pl)*KnEff/Kn
    P[init:,4] = np.random.normal(size=N-init)*np.sqrt(Pr/pr)*KnEff/Kn
    P[:init,5] = np.random.normal(size=init)*np.sqrt(Pl/pl)*KnEff/Kn
    P[init:,5] = np.random.normal(size=N-init)*np.sqrt(Pr/pr)*KnEff/Kn
    '''
    Tsim = 0
    return P, Tsim, N, Ne

#bins = 2*boxes
#hist = np.histogram(P[:,2],bins=bins,range=(0,1))[0]/N*bins*M 
#plt.bar(np.arange(bins)/bins+0.5/bins, hist, width = 1/boxes)
#plt.plot(np.arange(bins)/bins+0.5/bins, hist)
#plt.show()

def RunSim(sample, parallel=False):
    P, Tsim, N, Ne = Init(simcell)
    Tf = 0.1
    start = time.time()
    while Tsim < Tf:
        dT = np.min([dt, Tf-Tsim])
        print("Sample = ", sample,", Sim Time = ", Tsim, ", dt = ", dT, ", Real Time = ", HMS(time.time()-start) )

        #Advect
        ta = time.time()
        P = Advect(P,dT,dX) #This is dt not dteff
        print("Advected in in ", time.time()-ta)

        #Sort
        Boxes = Sort(P,boxes)

        #Collide
        tc = time.time()
        if parallel:
                Boxes = ParallelCollide(Boxes, dX, dT, s, Ne)
        else:
                Boxes = Collide(Boxes,dX,dT,s,Ne)
        print("Collided in ", time.time()-tc)

        P = np.concatenate(Boxes, axis=0)

        Tsim += dT

    bins = 100
    hist = np.histogram(P[:,2],bins=bins,range=(0,1))[0]
    rho = hist/N*bins*M
    vz = np.histogram(P[:,2], bins=bins, weights = P[:,5], range=(0,1))[0]/N*bins

    samplestr = str(sample)
    while len(samplestr) < 8:
        samplestr = '0' + samplestr
    np.save("SodData2/rho"+samplestr,rho)
    np.save("SodData2/vz"+samplestr,vz)
    end = time.time()
    print("Sample = ", sample, "Real Time = ", HMS(end-start) )

def ParallelSim(samples):
    if type(samples) == int:
	    samples = range(samples)
    elif type(samples) != type(range(1)):
        return("ParallelSim: Samples is not a list")
    p = mp.Pool()
    boxes = p.map(RunSim, samples)
    p.close()
    p.join()
    #return boxes


import sys

samples = range(int(sys.argv[1]), int(sys.argv[2]))

print(samples)
#ParallelSim(samples)
for sample in samples:
	RunSim(sample,True)
