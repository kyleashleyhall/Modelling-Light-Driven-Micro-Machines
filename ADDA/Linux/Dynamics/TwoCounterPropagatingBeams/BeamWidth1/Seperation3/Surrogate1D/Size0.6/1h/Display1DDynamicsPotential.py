# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants

PositionsFile = np.transpose(np.loadtxt('ParticlePositions'))

Positions=PositionsFile[1]

NumberOfBins=100

NumberOfIntermediatePoints=11 #Must be odd

plt.figure(1)
BinHeights,BinPositions=plt.hist(Positions,bins=NumberOfBins)[:2]
plt.show()
TotalArea=0

for i in range(NumberOfBins):
    TotalArea+=((BinPositions[i+1]-BinPositions[i])*BinHeights[i])
    
BinHeights/=TotalArea

kbT=(constants.Boltzmann)*((constants.zero_Celsius)+20)

PotentialPositions=np.zeros(0)

PotentialValues=np.zeros(0)

for i in range(0,NumberOfBins-2):
    
    ProbabilityStep=np.linspace(BinHeights[i],BinHeights[i+1],num=(NumberOfIntermediatePoints-1),endpoint=False)
    PotentialPositions=np.append(PotentialPositions,(np.linspace(BinPositions[i],BinPositions[i+1],num=(NumberOfIntermediatePoints-1),endpoint=False)))
    PotentialValues=np.append(PotentialValues,(0-(kbT*np.log(ProbabilityStep))))
    
ProbabilityStep=np.linspace(BinHeights[NumberOfBins-2],BinHeights[NumberOfBins-1],num=(NumberOfIntermediatePoints),endpoint=True)
PotentialPositions=np.append(PotentialPositions,(np.linspace(BinPositions[NumberOfBins-2],BinPositions[NumberOfBins-1],num=(NumberOfIntermediatePoints),endpoint=True)))
PotentialValues=np.append(PotentialValues,(0-(kbT*np.log(ProbabilityStep))))

PotentialPositions+=BinPositions[1]-BinPositions[0]

plt.figure(2)
plt.plot(PotentialPositions,PotentialValues)
plt.show()