# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:50:50 2018

@author: Kyle
"""

import numpy as np

def CalculateCentreOfMass(xPositions,yPositions,zPositions):
    
    x_CoM=(np.sum(xPositions))/(np.size(xPositions))
    y_CoM=(np.sum(yPositions))/(np.size(yPositions))
    z_CoM=(np.sum(zPositions))/(np.size(zPositions))
    
    CentreOfMass=np.array([[x_CoM],[y_CoM],[z_CoM]])
    
    return(CentreOfMass)
    
def CalculateMomentOfInertia(xPositions,yPositions,zPositions,CentreOfMass,Mass):
    
    DipMass=Mass/(np.size(xPositions))
    Inertia=np.zeros([3,1])
    Inertia[0,0]=np.sum((DipMass*((np.square(((yPositions-CentreOfMass[1,0])*(10**(-6)))))+(np.square(((zPositions-CentreOfMass[2,0])*(10**(-6))))))))
    Inertia[1,0]=np.sum((DipMass*((np.square(((xPositions-CentreOfMass[0,0])*(10**(-6)))))+(np.square(((zPositions-CentreOfMass[2,0])*(10**(-6))))))))
    Inertia[2,0]=np.sum((DipMass*((np.square(((xPositions-CentreOfMass[0,0])*(10**(-6)))))+(np.square(((yPositions-CentreOfMass[1,0])*(10**(-6))))))))
    
    return(Inertia)
    
def CalculateTorque(xPositions,yPositions,zPositions,F_x,F_y,F_z): #Assumes Forces are in Netwons!

    Torque=np.zeros([np.size(xPositions),6]) #Create an empty array to store position and the Torque in Nm at that point
    
    for i in range(np.size(xPositions)):
        Torque[i,0]=xPositions[i]
        Torque[i,1]=yPositions[i]
        Torque[i,2]=zPositions[i]
        Torque[i,3]=(((yPositions[i])*(F_z[i]))-((zPositions[i])*(F_y[i])))
        Torque[i,4]=(((zPositions[i])*(F_x[i]))-((xPositions[i])*(F_z[i])))
        Torque[i,5]=(((xPositions[i])*(F_y[i]))-((yPositions[i])*(F_x[i])))
        
    return(Torque)


CalculatedForces=np.transpose(np.loadtxt('CalculatedForces', skiprows=1))
Mass=1

CentreOfMass=CalculateCentreOfMass(IntFieldRaw[0,:],IntFieldRaw[1,:],IntFieldRaw[2,:])
Inertia=MomentOfInertia(IntFieldRaw[0,:],IntFieldRaw[1,:],IntFieldRaw[2,:],CentreOfMass,Mass)

print(CentreOfMass(IntFieldRaw[0,:],IntFieldRaw[1,:],IntFieldRaw[2,:]))