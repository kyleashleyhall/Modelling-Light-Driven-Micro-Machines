# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



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
    
    return np.vstack([axis,Force])
    
    
DipPolar = np.loadtxt('DipPol-X', skiprows=1) #skip first row if it is a heading
x, y, z = DipPolar[:,0], DipPolar[:,1], DipPolar[:,2]
DipPolx = DipPolar[:,4] + 1j*DipPolar[:,5]
DipPoly = DipPolar[:,6] + 1j*DipPolar[:,7]
DipPolz = DipPolar[:,8] + 1j*DipPolar[:,9]

Eintn = np.loadtxt('IntField-X', skiprows=1)
Eintx = Eintn[:,4] + 1j*Eintn[:,5]
Einty = Eintn[:,6] + 1j*Eintn[:,7]
Eintz = Eintn[:,8] + 1j*Eintn[:,9]

axis = np.vstack([x,y,z])
entry = np.vstack([DipPolx, DipPoly, DipPolz])
yentr = np.vstack([Eintx, Einty, Eintz])

#print(Forces(axis, entry, yentr))
Quiv = Forces(axis, entry, yentr)

plt.quiver(Quiv[0], Quiv[1], Quiv[4], Quiv[5])
plt.show()


print(sum(Quiv[3]), sum(Quiv[4]), sum(Quiv[5]))











