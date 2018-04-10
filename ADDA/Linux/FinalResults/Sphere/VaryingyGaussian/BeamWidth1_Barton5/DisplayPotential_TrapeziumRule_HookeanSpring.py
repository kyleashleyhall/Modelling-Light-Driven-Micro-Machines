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
            
beamwidth=1
impedance= (376.73/1.57)
t0 = np.linspace(-5,5,len(Positions))
minval = np.amin(TrapeziumPotential)
cons = 0.5
cons_step = 0.1
r = minval*np.exp(-cons*(t0**2))[35:65]
TrapPotential = TrapeziumPotential[35:65]
costval = np.sum((TrapPotential-r)**2)
print(costval)
udtracker = np.array([1])
while costval > 8e-29:
#for i in range(3):
    consU = cons + cons_step
    rU = minval*np.exp(-consU*(t0**2))[35:65]
    costvalU = np.sum((TrapPotential-rU)**2)
    consD = cons - cons_step
    rD = minval*np.exp(-consD*(t0**2))[35:65]
    costvalD = np.sum((TrapPotential-rD)**2)
    if costvalU < costvalD:
        print('yes')
        cons += cons_step
        costval = costvalU
        udtracker = np.append(udtracker, np.array([1]), axis=0)
    else:
        print('no')
        cons -= cons_step
        costval = costvalD
        udtracker = np.append(udtracker, np.array([0]), axis=0)
        
    if np.sum(udtracker[-2:]) == 1:
        cons_step *= 0.5
print(costval)
print(cons)
    
Analtical = minval*np.exp(-cons*(t0**2))


plt.figure(1)

TrapeziumPlot=plt.errorbar(Positions,TrapeziumPotential,label="Trapezium Rule Potential")
HookeanPlot=plt.errorbar(HookeanPositions,HookeanPotential,label="Hookean Potential")
plt.plot(t0,Analtical, label='2')
plt.legend(handles=[TrapeziumPlot,HookeanPlot])

plt.show()