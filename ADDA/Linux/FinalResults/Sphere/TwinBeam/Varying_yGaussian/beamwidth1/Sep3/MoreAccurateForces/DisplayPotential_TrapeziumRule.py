# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Forces = np.transpose(np.loadtxt('TotalForce'))

Potential=np.zeros(np.size(Forces,axis=1))

Positions=Forces[0]

Force=Forces[1]

Seperation=Positions[1]-Positions[0]

Potential[1]=0-((Seperation/2)*Force[1])

for iterator1 in range(2,np.size(Forces,axis=1)):
    
    for iterator2 in range(1,iterator1):
        
        Potential[iterator1]-=(Seperation*Force[iterator2])
        
    Potential[iterator1]-=((Seperation/2)*Force[iterator1])
        
    

plt.figure(1)

plt.plot(Positions,Potential)

plt.show()