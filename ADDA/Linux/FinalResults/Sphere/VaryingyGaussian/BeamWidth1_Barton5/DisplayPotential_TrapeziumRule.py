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
        
beamwidth=1
impedance= (376.73/1.57)
t0 = np.linspace(-5,5,10000)
minval = np.amin(Potential)
print(minval)
cons = minval*np.pi*((2*beamwidth*10**-6)**2)/(2*impedance*5e-3)
print(cons)
r = minval*np.exp(-(t0**2)/(beamwidth**2))



plt.figure(1)

plt.plot(Positions,Potential, label='trap')
plt.plot(t0,r, label='2')
plt.legend()
plt.show()