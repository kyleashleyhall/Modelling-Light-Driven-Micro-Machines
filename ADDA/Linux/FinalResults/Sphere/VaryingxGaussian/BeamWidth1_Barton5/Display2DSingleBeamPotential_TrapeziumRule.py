# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

xForcesFile = np.transpose(np.loadtxt('..\..\VaryingxGaussian\BeamWidth1_Barton5\TotalForce'))
yForcesFile = np.transpose(np.loadtxt('..\..\VaryingyGaussian\BeamWidth1_Barton5\TotalForce'))

Potential=np.zeros([np.size(yForcesFile,axis=1),np.size(xForcesFile,axis=1)])
xTempPotential=np.zeros(np.size(xForcesFile,axis=1))
yTempPotential=np.zeros(np.size(yForcesFile,axis=1))

xPositions=xForcesFile[0]
xForce=xForcesFile[3]

yPositions=yForcesFile[1]
yForce=yForcesFile[4]

xSeperation=xPositions[1]-xPositions[0]
ySeperation=xPositions[1]-xPositions[0]

xZeroArgument=np.argwhere(np.isclose(xPositions,0))[0,0]
yZeroArgument=np.argwhere(np.isclose(yPositions,0))[0,0]

xTempPotential[1]=0-((xSeperation/2)*xForce[1])

for iterator1 in range(2,np.size(xForce)):
    
    for iterator2 in range(1,iterator1):
        
        xTempPotential[iterator1]-=(xSeperation*xForce[iterator2])
        
    xTempPotential[iterator1]-=((xSeperation/2)*xForce[iterator1])
    
yTempPotential[1]=0-((ySeperation/2)*yForce[1])

for iterator1 in range(2,np.size(yForce)):
    
    for iterator2 in range(1,iterator1):
        
        yTempPotential[iterator1]-=(ySeperation*yForce[iterator2])
        
    yTempPotential[iterator1]-=((ySeperation/2)*yForce[iterator1])

Potential[:,xZeroArgument]=yTempPotential
Potential[yZeroArgument,:]=xTempPotential

for j in range(np.size(Potential,axis=0)):
    for i in range(np.size(Potential,axis=1)):
        
        if(i!=xZeroArgument):
            if(j!=yZeroArgument):
                R_value=np.sqrt(((np.square(xPositions[i]))+(np.square(yPositions[j]))))
                if((np.any(xPositions>R_value))&(np.any(yPositions>R_value))):
                    xArgumentAboveR=np.argwhere(xPositions>R_value)[0,0]
                    yArgumentAboveR=np.argwhere(yPositions>R_value)[0,0]
                    xPotentialBelowR=Potential[xArgumentAboveR-1,0]
                    yPotentialBelowR=Potential[0,yArgumentAboveR-1]
                    xPotentialAtR=xTempPotential[xArgumentAboveR-1]+((xTempPotential[xArgumentAboveR]-xTempPotential[xArgumentAboveR-1])*((R_value-xPositions[xArgumentAboveR-1])/(xPositions[xArgumentAboveR]-xPositions[xArgumentAboveR-1])))
                    yPotentialAtR=yTempPotential[yArgumentAboveR-1]+((yTempPotential[yArgumentAboveR]-yTempPotential[yArgumentAboveR-1])*((R_value-yPositions[yArgumentAboveR-1])/(yPositions[yArgumentAboveR]-yPositions[yArgumentAboveR-1])))
                    Theta=np.arccos(xPositions[i]/R_value)
                    Potential[j,i]=0.5*((yPotentialAtR+xPotentialAtR)+((yPotentialAtR-xPotentialAtR)*(np.sin((2*Theta)-(0.5*np.pi)))))
                
                else:
                    Potential[j,i]=0
                
        
        

xPositions-=(xPositions[1]-xPositions[0])/2
yPositions-=(yPositions[1]-yPositions[0])/2

plt.figure(1)
plt.pcolormesh(xPositions,yPositions,np.transpose(Potential),cmap='bone')
plt.colorbar()
plt.show()