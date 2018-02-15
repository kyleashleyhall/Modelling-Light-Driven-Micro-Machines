# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:01:21 2018

@author: Kyle
"""

import numpy as np

def CalculateMomentOfInertia(xPositions,yPositions,zPositions,CentreOfMass,Mass):
    
    DipMass=Mass/(np.size(xPositions))
    Inertia=np.zeros([3,1])
    Inertia[0,0]=np.sum((DipMass*((np.square(((yPositions-CentreOfMass[1,0])*(10**(-6)))))+(np.square(((zPositions-CentreOfMass[2,0])*(10**(-6))))))))
    Inertia[1,0]=np.sum((DipMass*((np.square(((xPositions-CentreOfMass[0,0])*(10**(-6)))))+(np.square(((zPositions-CentreOfMass[2,0])*(10**(-6))))))))
    Inertia[2,0]=np.sum((DipMass*((np.square(((xPositions-CentreOfMass[0,0])*(10**(-6)))))+(np.square(((yPositions-CentreOfMass[1,0])*(10**(-6))))))))
    
    return(Inertia)
    
IntFieldRaw=np.transpose(np.loadtxt('IntField-Y', skiprows=1))
CentreOfMass=np.array([[0],[0],[0]])
Mass=1
Inertia=MomentOfInertia(IntFieldRaw[0,:],IntFieldRaw[1,:],IntFieldRaw[2,:],CentreOfMass,Mass)

print(Inertia[0,0])
print(Inertia[1,0])
print(Inertia[2,0])