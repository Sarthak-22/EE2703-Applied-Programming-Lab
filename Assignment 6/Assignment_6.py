 
# Sarthak Vora - EE19B140 (Asignment 6)

# Run the code using the following command - 'python EE19B140_assgn6.py n=value M=value Msig=Value nk=Value u0=Value p=Value seed=Value'

import numpy as np    
import matplotlib.pyplot as plt
import pandas as pd    # pandas for showing in tabular form
import sys


n = 100    # Spatial grid size
M = 5     # no. of electrons injected per turn
Msig = 2   # standard deviation of the distribution electrons injected	
nk = 500   # no. of turns to simulate
u0 = 5.0   # threshold velocity
p = 0.25   # probability that ionization will occur
seedValue = 2   # A seed for random number generator


args = sys.argv[1:]  # This is because argv[0] is the program name itself which is not required
argList = {}         # Dictionary declared
for ele in args:
    ele = ele.split('=')    # 'varname=value' string is split into varname and value
    argList[ele[0]] = float(ele[1])   # The dictionary is updated with values from commandline


if 'n' in argList:
    n = argList['n']
if 'M' in argList:
    M = argList['M']
if 'nk' in argList:
    nk = argList['nk']
if 'u0' in argList:
    u0 = argList['u0']
if 'p' in argList:
    p = argList['p']
if 'Msig' in argList:
    Msig = argList['Msig']
if 'seed' in argList:
    seedValue = argList['seed']


xx = np.zeros(int(n*M))  # electron position
u = np.zeros(int(n*M))   # electron velocity
dx = np.zeros(int(n*M))  # electron displacement
np.random.seed((int)(seedValue))   # A seed given so that uniform numbers results across runs for the random values generated in the code

I = []  # This is used to store the photons generated at each location in each turn
X = []  # This stores the position of every electron after each turn
V = []  # Stores the velocity of the electron after each turn


def electronCollision(u,xx,dx,kl):    # This function does the collision update as per the question    
    u[kl] = 0.0
    xx[kl] = xx[kl] - dx[kl]*np.random.rand(1)[0]


def electronCollisionModified(u,xx,dx,kl):   # This function does the collision update more accurately taking a random distribution of time 
    t = np.random.rand(1)[0]
    xx[kl] = xx[kl] - dx[kl]
    u[kl] = u[kl] - 1.0
    xx[kl] = xx[kl] + u[kl]*t + 0.5*t*t
    u[kl] = 0.0


ii = []
for k in range(1,int(nk)):
    dx[ii] = u[ii] + 0.5 
    xx[ii] = xx[ii] + dx[ii]
    u[ii] = u[ii] + 1.0
    
    jj = np.where(xx > n)[0]
    dx[jj] = 0.0
    xx[jj] = 0.0
    u[jj] = 0.0
    
    kk = np.where( u >= u0 )[0]
    ll = np.where(np.random.rand(len(kk)) <= p)
    kl = kk[ll]

    electronCollisionModified(u, xx, dx, kl)
    I.extend(xx[kl].tolist())
    
    m = np.random.randn()*Msig + M     
    ll = np.where(xx == 0)[0]
    maxElec = min(len(ll),(int)(m))
    xx[ll[0:maxElec]] = 1.0
    u[ll[0:maxElec]] = 0.0

    ii = np.where(xx > 0)[0]
    X.extend(xx[ii].tolist())
    V.extend(u[ii].tolist())



# Population plot for electron density
plt.hist(X,histtype='bar',range=(10,n), bins=np.arange(1,n,n/100),ec='black',alpha=0.5) 
plt.title('Population Plot')
plt.xlabel('x')
plt.ylabel('No. of electrons')
plt.show()


# Population plot for intensity of emitted light
plt.hist(I,histtype='bar',range=(10,n), bins=np.arange(1,n,n/100),ec='black',alpha=0.5)
plt.title('Intensity Plot')
plt.xlabel('x')
plt.ylabel('No. of photons emitted ($\propto$ Intensity)')
plt.show()


# Electron Phase plot
plt.plot(X,V,'x')
plt.title('Electron phase space')
plt.xlabel('x')
plt.ylabel('velocity')
plt.show()


# Finding the intensity table
bins = plt.hist(I,bins=np.arange(1,n,n/100))[1]    # Bin positions are obtained
count = plt.hist(I,bins=np.arange(1,n,n/100))[0]   # Population counts obtained
xpos = 0.5*(bins[0:-1] + bins[1:])     # As no. of end-points of bins would be 1 more than actual no. of bins, the mean of bin end-points are used to get population of count a particular bin
df = pd.DataFrame()   # A pandas dataframe is initialized to do the tabular plotting of values.
df['Intensity'] = xpos
df['Count'] = count
print (df)