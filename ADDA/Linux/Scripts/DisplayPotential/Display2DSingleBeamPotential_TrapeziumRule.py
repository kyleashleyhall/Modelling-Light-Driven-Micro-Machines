# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

xForcesFile = np.transpose(np.loadtxt('TotalForce'))
yForcesFile = np.transpose(np.loadtxt('TotalForce'))

Potential=np.zeros([np.size(yForcesFile,axis=1),np.size(xForcesFile,axis=1)])
xTempPotential=np.zeros(np.size(xForcesFile,axis=1))
yTempPotential=np.zeros(np.size(yForcesFile,axis=1))

xPositions=xForcesFile[1]
xForce=xForcesFile[4]

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
    
#Set the minimum of the x and y equal to 0

xTempPotential-=np.amin(xTempPotential)
yTempPotential-=np.amin(yTempPotential)


Potential[:,xZeroArgument]=yTempPotential
Potential[yZeroArgument,:]=xTempPotential

for j in range(np.size(yForcesFile,axis=1)):
    for i in range(np.size(xForcesFile,axis=1)):
        
        while(j!=yZeroArgument):
            while(i!=xZeroArgument):
                Potential[j,i]=np.sqrt(np.power((yTempPotential[j]),2)+np.power((xTempPotential[i]),2))
        
xPositions-=(xPositions[1]-xPositions[0])/2
yPositions-=(yPositions[1]-yPositions[0])/2

plt.figure(1)
plt.pcolormesh(xPositions,yPositions,np.transpose(Potential),cmap='bone')
plt.colorbar()
plt.show()