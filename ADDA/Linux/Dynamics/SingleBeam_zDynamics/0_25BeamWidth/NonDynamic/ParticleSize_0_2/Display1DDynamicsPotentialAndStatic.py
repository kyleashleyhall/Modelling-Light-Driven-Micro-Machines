# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants
import os

ParticleSize1=int(os.getcwd()[-3])
ParticleSize2=int(os.getcwd()[-1])
DynamicPositionsFile = np.transpose(np.loadtxt('..'+os.sep+'..'+os.sep+'z_0Start'+os.sep+'ParticleSize_'+str(ParticleSize1)+'_'+str(ParticleSize2)+os.sep+'ParticlePositions'))

Forces = np.transpose(np.loadtxt('..'+os.sep+'..'+os.sep+'NonDynamic'+os.sep+'ParticleSize_'+str(ParticleSize1)+'_'+str(ParticleSize2)+os.sep+'TotalForce'))

Potential=np.zeros(np.size(Forces,axis=1))

StaticPositions=Forces[2]

Force=Forces[5]

Seperation=StaticPositions[1]-StaticPositions[0]

Potential[1]=0-((Seperation/2)*Force[1])

for iterator1 in range(2,np.size(Forces,axis=1)):
    
    for iterator2 in range(1,iterator1):
        
        Potential[iterator1]-=(Seperation*Force[iterator2])
        
    Potential[iterator1]-=((Seperation/2)*Force[iterator1])

Positions=DynamicPositionsFile[3]

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
DynamicPlot=plt.errorbar(PotentialPositions,PotentialValues,label="Dynamic Potential")
StaticPlot=plt.errorbar(StaticPositions,Potential,label="Statics Potential")
plt.legend(handles=[DynamicPlot,StaticPlot])
plt.show()