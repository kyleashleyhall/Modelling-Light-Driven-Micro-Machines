# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:41:45 2018

@author: Kyle
"""

import numpy as np

def CalculateCentreOfMass(xPositions,yPositions,zPositions):
    
    x_CoM=(np.sum(xPositions))/(np.size(xPositions))
    y_CoM=(np.sum(yPositions))/(np.size(yPositions))
    z_CoM=(np.sum(zPositions))/(np.size(zPositions))
    
    CentreOfMass=np.array([[x_CoM],[y_CoM],[z_CoM]])
    
    return(CentreOfMass)

IntFieldRaw=np.transpose(np.loadtxt('IntField-Y', skiprows=1))
print(CentreOfMass(IntFieldRaw[0,:],IntFieldRaw[1,:],IntFieldRaw[2,:]))