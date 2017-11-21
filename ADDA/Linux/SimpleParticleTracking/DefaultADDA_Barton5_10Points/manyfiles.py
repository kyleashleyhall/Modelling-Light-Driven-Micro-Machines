# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import glob
import os



def Forces(axes, DipPol, Eint):
    d_l = np.zeros([3])  #d for central difference formula, for each axis
    for i in range(len(d_l)):
        d_l[i] = (max(axes[i]) - min(axes[i])) / len(axes[i])
        
    dE = np.zeros([3, 3, len(Eint[0])], dtype=complex) #Derivative at each point, for each axis
    Force = np.zeros([3, len(Eint[0])])
    for i in range(3):
        d = d_l[i]
        for k in range(3):
            for j in range(len(Eint[i])-2):
                dE[i, k, j+1] = (Eint[k, j+2] - Eint[k, j]) / (2 * d) #Central difference formula
            
            dE[i,k,0], dE[i,k,-1] = dE[i,k,1], dE[i,k,-2] #set the boundary values
        
        Force[i] =  np.real(DipPol[i]*np.conjugate(sum(dE[i]))) /2 #Calculate the Force
    
    return np.vstack([axes,Force])
    
def FileSlice(Fname):
    filename = np.loadtxt(Fname, skiprows=1)    
    x,y,z = filename[:,0], filename[:,1], filename[:,2]
    varX = filename[:,4] + 1j*filename[:,5]
    varY = filename[:,6] + 1j*filename[:,7]
    varZ = filename[:,8] + 1j*filename[:,9]
    return np.vstack([x,y,z]), np.vstack([varX,varY,varZ])


print('---------------------')
DipPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'DipPol-X')
IntFPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IntField-X')
DipFiles, IntFFiles = glob.glob(DipPathInput), glob.glob(IntFPathInput) #File containing the paths to each DipPol, IntField file
for i in range(len(DipFiles)): 
    axes = FileSlice(DipFiles[i])[0] #array containing the position of all dipoles
    DipPol = FileSlice(DipFiles[i])[1] #array with the polarization of said dipoles
    IntField = FileSlice(IntFFiles[i])[1] #array with the Internal field at each point
    Force = Forces(axes, DipPol, IntField) #Calculation of the force at each point
    #Add analysis here eg quiver plot, total force calculation etc to be done on each simulation
    plt.quiver(Force[0], Force[1], Force[2], Force[3])
    plt.show()
    
#Add any overall analysis here










