# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
import os

Forces = np.transpose(np.loadtxt('TotalForce'))

GridSize=int(np.sqrt(np.size(Forces,axis=1)))

xPositions=np.zeros(GridSize)

yPositions=np.zeros(GridSize)

xForceGrid=np.zeros([GridSize,GridSize])

yForceGrid=np.zeros([GridSize,GridSize])

Potential=np.zeros([GridSize,GridSize])

for i in range (GridSize):
    for j in range (GridSize):
        iterator=(i*GridSize)+j
        if(j==0):
            xPositions[i]=Forces[0,iterator]
        if(i==0):
            yPositions[j]=Forces[1,iterator]
            
        xForceGrid[j,i]=Forces[3,iterator]
        yForceGrid[j,i]=Forces[4,iterator]

xSeperation=xPositions[1]-xPositions[0]
ySeperation=yPositions[1]-yPositions[0]

#Initial y Integral

Potential[1,0]=0-((ySeperation/2)*yForceGrid[1,0])

for iterator1 in range(2,GridSize):
    
    for iterator2 in range(1,iterator1):
        
        Potential[iterator1,0]-=(ySeperation*yForceGrid[iterator2,0])
        
    Potential[iterator1,0]-=((ySeperation/2)*yForceGrid[iterator1,0])
    
#Initial x Integral
    
for iterator0 in range(GridSize):
    
    Potential[iterator0,1]=Potential[iterator0,0]-((xSeperation/2)*xForceGrid[iterator0,1])
    
    for iterator1 in range(2,GridSize):
    
        for iterator2 in range(1,iterator1):
        
            Potential[iterator0,iterator1]-=(xSeperation*xForceGrid[iterator0,iterator2])
        
        Potential[iterator0,iterator1]-=((xSeperation/2)*xForceGrid[iterator0,iterator2])
        #Potential[iterator0,iterator1]+=Potential[iterator0,0]
        
BeamWidth=int(((os.getcwd())[-4:-3]))+((int(((os.getcwd())[-2:]))/100))
TitleString=r"Torque about the z axis on a cone in the xy plane for a "+str(BeamWidth)+" $\mu m$ beam width"

plt.figure(1)
ax = plt.gca()
ax.set_aspect('equal', 'box')
plt.title(TitleString)
plt.xlabel(r"Particle position, x ($\mu m$)")
plt.ylabel(r"Particle position, y ($\mu m$)")
PotentialPlot=plt.pcolormesh(xPositions,yPositions,Potential,cmap='bone',vmin=-1e-20,vmax=1e-20)
BeamCircle=plt.Circle((0, 0), BeamWidth, color='b', fill=False,label='Beam Focus')
ParticleInBeamCircle=plt.Circle((0, 0), (BeamWidth+1), color='r', fill=False,label='Particle in Beam')
ax.add_artist(BeamCircle)
ax.add_artist(ParticleInBeamCircle)
plt.legend(handles=[BeamCircle,ParticleInBeamCircle])
Colorbar=plt.colorbar()
Colorbar.set_label("Potential")
plt.show()