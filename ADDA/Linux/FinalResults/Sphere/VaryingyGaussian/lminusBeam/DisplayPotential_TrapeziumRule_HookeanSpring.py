# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Forces = np.transpose(np.loadtxt('TotalForce'))

TrapeziumPotential=np.zeros(np.size(Forces,axis=1))

Positions=Forces[1]

Force=Forces[4]

Seperation=Positions[1]-Positions[0]

TrapeziumPotential[1]=0-((Seperation/2)*Force[1])

for iterator1 in range(2,np.size(Forces,axis=1)):
    
    for iterator2 in range(1,iterator1):
        
        TrapeziumPotential[iterator1]-=(Seperation*Force[iterator2])
        
    TrapeziumPotential[iterator1]-=((Seperation/2)*Force[iterator1])
        
MinPotentialArgument=np.argmin(TrapeziumPotential)

kValue=(Force[MinPotentialArgument-1]-Force[MinPotentialArgument+1])/(2*Seperation)
print(kValue)

TempHookeanPotential=(0.5*kValue*(np.power((Positions-Positions[MinPotentialArgument]),2)))+np.amin(TrapeziumPotential)

HookeanPotential=np.zeros(0)
HookeanPositions=np.zeros(0)

MaxPotentialValue=0
PositionIdentifier=0
for iterator1 in range(np.size(Forces,axis=1)):
    
    if(TempHookeanPotential[iterator1]<MaxPotentialValue):
        if(PositionIdentifier==0):
            HookeanPotential=np.append(HookeanPotential,TempHookeanPotential[iterator1-1])
            HookeanPositions=np.append(HookeanPositions,Positions[iterator1-1])
            PositionIdentifier=1
        
        HookeanPotential=np.append(HookeanPotential,TempHookeanPotential[iterator1])
        HookeanPositions=np.append(HookeanPositions,Positions[iterator1])
        
    else:
        if(PositionIdentifier==1):
            HookeanPotential=np.append(HookeanPotential,TempHookeanPotential[iterator1])
            HookeanPositions=np.append(HookeanPositions,Positions[iterator1])
            PositionIdentifier=0

plt.figure(1)

TrapeziumPlot=plt.errorbar(Positions,TrapeziumPotential,label="Trapezium Rule Potential")
HookeanPlot=plt.errorbar(HookeanPositions,HookeanPotential,label="Hookean Potential")

plt.legend(handles=[TrapeziumPlot,HookeanPlot])

plt.show()