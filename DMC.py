from scipy.interpolate import griddata
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np

#Shell environment
import os
import sys
import subprocess		# to send python variables to shell.

# Permutations
import itertools

# He

def distance(_x):
    return np.sqrt(np.dot(_x,_x))

#Change here.
def potentialV(_xp):                   
#     He atom                        
#    r1 = distance(_xp[0:nD])
#    r2 = distance(_xp[nD:2*nD])
#    r12 = distance(_xp[0:nD] - _xp[nD:2*nD])
#    return -2.0/r1 - 2.0/r2 + 1.0/r12
#     H^- ion                        
    r1 = distance(_xp[0:nD])
    r2 = distance(_xp[nD:2*nD])
    r12 = distance(_xp[0:nD] - _xp[nD:2*nD])
    return -1.0/r1 - 1.0/r2 + 1.0/r12
#     H- atom                        
#    r1 = distance(_xp[0:nD])
#    return -1.0/r1 

### Initiate variables
nP = 2 # Number of particles
nD = 3 # Number of dimension
nM = 1000 # Number of walkers.
Alpha = 1.0
dtau = 0.02
nDMC = 1000
###

xsample = np.array( [ np.random.uniform(-10.0,10.0, nD*nP) for _iM in range(nM)]) # Initial random starting point.
# xsample = np.array( [ np.linspace(-5.0,5.0, nD*nP) for _iM in range(nM)]) # Initial random starting point.

#list_Eguess = np.array([]).reshape(0,4)
f_para = open('./dmcout.dat','w')

for _iDMC in range(nDMC):
    _xproposed = np.array([]).reshape(0,nP*nD)
    _nM = len(xsample)
    
    Eguess = 0.0
    for _iP in range(_nM):
        Eguess = Eguess + potentialV(xsample[_iP])/_nM
    Eguess = Eguess - Alpha*np.log(_nM/nM)
    #list_Eguess = np.vstack(( list_Eguess, np.array([_iDMC, Eguess, _nM, nM]) ))
    np.savetxt(f_para, np.array([[_iDMC, Eguess, _nM, nM]]) , delimiter=" ", fmt="%0.12f")

    xsample = xsample +  np.array( [np.random.normal(0.0, np.sqrt(dtau) ,nD*nP) for _iM in range(_nM)] )
    
    for _iM in range(_nM):
         # Branching.
        _potential = np.exp( -dtau*( potentialV(xsample[_iM]) - Eguess ) ) 
        _mNumber = int( _potential + np.random.uniform(0,1,1) )
        if _mNumber > 0 :
            _xproposed = np.append(_xproposed, [xsample[_iM]], axis=0)  
        if _mNumber > 1 :
            _xproposed = np.append(_xproposed, [xsample[_iM]], axis=0)
    xsample = np.copy(_xproposed)

#####
