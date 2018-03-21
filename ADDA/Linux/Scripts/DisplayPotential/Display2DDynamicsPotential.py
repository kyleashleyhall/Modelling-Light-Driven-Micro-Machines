# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants

PositionsFile = np.transpose(np.loadtxt('ParticlePositions'))

Positions1=PositionsFile[1]
Positions2=PositionsFile[2]

NumberOfBins=100

NumberOfIntermediatePoints=11 #Must be odd

plt.figure(1)
BinHeights,xBinPositions,yBinPositions=plt.hist2d(Positions1,Positions2,bins=NumberOfBins,cmap='bone')[:3]
plt.colorbar()
plt.show()
TotalVolume=0

for i in range(NumberOfBins):
    for j in range(NumberOfBins):
        TotalVolume+=((xBinPositions[i+1]-xBinPositions[i])*(yBinPositions[j+1]-yBinPositions[j])*BinHeights[i,j])
        
print(BinHeights)
BinHeights/=TotalVolume

kbT=(constants.Boltzmann)*((constants.zero_Celsius)+20)

PotentialxPositions=(xBinPositions)
PotentialyPositions=(yBinPositions)
PotentialGrid=0-(kbT*np.log(BinHeights))

for i in range(NumberOfBins):
    for j in range(NumberOfBins):
        if(PotentialGrid[i,j]==np.infty):
            PotentialGrid[i,j]=-np.infty
            
MaxPotential=np.amax(PotentialGrid)

for i in range(NumberOfBins):
    for j in range(NumberOfBins):
        if(PotentialGrid[i,j]==-np.infty):
            PotentialGrid[i,j]=MaxPotential

plt.figure(2)
plt.pcolormesh(PotentialxPositions,PotentialyPositions,np.transpose(PotentialGrid),cmap='bone')
plt.colorbar()
plt.show()