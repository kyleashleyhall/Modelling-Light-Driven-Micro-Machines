# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Forces = np.transpose(np.loadtxt('TotalForce'))

Potential=np.zeros(np.size(Forces,axis=1))

Positions=Forces[1]

Force=Forces[4]

Seperation=Positions[1]-Positions[0]

Potential[1]=0-((Seperation/2)*Force[1])

for iterator1 in range(2,np.size(Forces,axis=1)):
    
    for iterator2 in range(1,iterator1):
        
        Potential[iterator1]-=(Seperation*Force[iterator2])
        
    Potential[iterator1]-=((Seperation/2)*Force[iterator1])
        
Potential-=np.amin(Potential)

PotentialZeroArgument=np.argmin(Potential)

Variable_k=np.zeros(np.size(Potential))

for iterator1 in range(np.size(Variable_k)):
    
    if(iterator1!=PotentialZeroArgument):
        Variable_k[iterator1]=(2*Potential[iterator1])/(np.power((Positions[iterator1]-Positions[PotentialZeroArgument]),2))
        
Variable_k[PotentialZeroArgument]=(Variable_k[(PotentialZeroArgument+1)]+Variable_k[(PotentialZeroArgument-1)])/2

plt.figure(1)

TrapPlot=plt.errorbar(Positions,Potential,label="Trapezium Rule")
kPlot=plt.errorbar(Positions,(0.5*np.multiply(Variable_k,(np.square(Positions-Positions[PotentialZeroArgument])))),label=r"Variable $k(x)$ interpretation")

plt.legend(handles=[TrapPlot,kPlot])

plt.show()