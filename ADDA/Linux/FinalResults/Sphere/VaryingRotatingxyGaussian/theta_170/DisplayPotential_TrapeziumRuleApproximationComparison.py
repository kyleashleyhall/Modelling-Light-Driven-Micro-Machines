# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

theta=int(os.getcwd()[-3:])

Forces = np.transpose(np.loadtxt('TotalForce'))

Potential=np.zeros(np.size(Forces,axis=1))

Positions=Forces[3]

Force=Forces[7]

Seperation=Positions[1]-Positions[0]

Potential[1]=0-((Seperation/2)*Force[1])

for iterator1 in range(2,np.size(Forces,axis=1)):
    
    for iterator2 in range(1,iterator1):
        
        Potential[iterator1]-=(Seperation*Force[iterator2])
        
    Potential[iterator1]-=((Seperation/2)*Force[iterator1])
        
xForcesFile=np.transpose(np.loadtxt('..\theta_000\TotalForce'))
xPotential=np.zeros(np.size(xForcesFile,axis=1))
xPositions=xForcesFile[3]
xForce=xForcesFile[7]
xSeperation=xPositions[1]-xPositions[0]

xPotential[1]=0-((xSeperation/2)*xForce[1])

for iterator1 in range(2,np.size(xForcesFile,axis=1)):
    
    for iterator2 in range(1,iterator1):
        
        xPotential[iterator1]-=(xSeperation*xForce[iterator2])
        
    xPotential[iterator1]-=((xSeperation/2)*xForce[iterator1])

yForcesFile=np.transpose(np.loadtxt('..\theta_090\TotalForce'))
yPotential=np.zeros(np.size(yForcesFile,axis=1))
yPositions=yForcesFile[3]
yForce=yForcesFile[7]
ySeperation=yPositions[1]-yPositions[0]

yPotential[1]=0-((ySeperation/2)*yForce[1])

for iterator1 in range(2,np.size(yForcesFile,axis=1)):
    
    for iterator2 in range(1,iterator1):
        
        yPotential[iterator1]-=(ySeperation*yForce[iterator2])
        
    yPotential[iterator1]-=((ySeperation/2)*yForce[iterator1])
    
ApproximationPotential=np.zeros(np.size(Potential))

for iterator1 in range(np.size(Potential)):
    
    xValue=Positions[iterator1]*np.cos(theta)
    yValue=Positions[iterator1]*np.sin(theta)
    if(np.any(xPositions>xValue)):
        xArgumentAboveR=np.argwhere(xPositions>xValue)[0,0]
        if(xArgumentAboveR==0):
            xPotential[iterator1]=0
        else:
            xPotential[iterator1]=xPotential[xArgumentAboveR-1]+((xPotential[xArgumentAboveR]-xPotential[xArgumentAboveR-1])*((xValue-xPositions[xArgumentAboveR-1])/(xPositions[xArgumentAboveR]-xPositions[xArgumentAboveR-1])))
    else:
        xPotential[iterator1]=0
        
    if(np.any(yPositions>yValue)):
        yArgumentAboveR=np.argwhere(yPositions>yValue)[0,0]
        if(yArgumentAboveR==0):
            yPotential[iterator1]=0
        else:
            yPotential[iterator1]=yPotential[yArgumentAboveR-1]+((yPotential[yArgumentAboveR]-yPotential[yArgumentAboveR-1])*((yValue-yPositions[yArgumentAboveR-1])/(yPositions[yArgumentAboveR]-yPositions[yArgumentAboveR-1])))
    else:
        yPotential[iterator1]=0
        
    ApproximationPotential[iterator1]=(0.5*(xPotential[iterator1]+xPotential[iterator1]))+((yPotential[iterator1]-xPotential[iterator1])*(np.sin((2*theta)-(0.5*np.pi))))

plt.figure(1)

AccuratePlot=plt.errorbar(Positions,Potential,label="Trapezium Potential")
ApproximationPlot=plt.errorbar(Positions,Potential,label="Approximation")
plt.legend(handles=[AccuratePlot,ApproximationPlot])

plt.show()