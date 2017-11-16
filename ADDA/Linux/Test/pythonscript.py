# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt

def RadForce(x, y, z, l, DipPoll, Eintl):  #l is an axis, so either x, y or z.
    d = (max(l) - min(l)) / len(l)  #takes the width of the shape, and divides it by the number of dipoles CHECK
    dE = np.empty([len(Eintl),len(Eintl)], dtype=complex) 
    for i in range(len(Eintl)-2):
        dE[i+1] = ((Eintl[i] - Eintl[i+2]))/ (2 * d)  #Central Difference Formula with error in d^2 to calculate the derivative
    dE[0] = dE[1] #set first value due to boundary
    dE[-1] = dE[-2] #set last value due to boundary

    Force = 0.5 * np.real(DipPoll*np.conjugate(dE))
    
    return np.vstack([x,y,z,Force])

    
DipPolar = np.loadtxt('DipPol-X', skiprows=1) #skip first row if it is a heading
x, y, z = DipPolar[:,0], DipPolar[:,1], DipPolar[:,2]
DipPolx = DipPolar[:,4] + 1j*DipPolar[:,5]
DipPoly = DipPolar[:,6] + 1j*DipPolar[:,7]
DipPolz = DipPolar[:,8] + 1j*DipPolar[:,9]

Eintn = np.loadtxt('IntField-X', skiprows=1)
Eintx = Eintn[:,4] + 1j*Eintn[:,5]
Einty = Eintn[:,6] + 1j*Eintn[:,7]
Eintz = Eintn[:,8] + 1j*Eintn[:,9]


print(RadForce(x,y,z,y,DipPolx,Eintx))




#Going to need all 3 vectors for output










