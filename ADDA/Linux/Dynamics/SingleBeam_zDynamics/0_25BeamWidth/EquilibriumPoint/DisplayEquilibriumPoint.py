# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
import os

ParticleSize=np.arange(0.2,2.2,0.2)
EquilibriumPoint=np.zeros(np.size(ParticleSize))

for iterator in range(np.size(ParticleSize)):
    LastDigit=((iterator*2)+2)%10
    FirstDigit=((iterator*2)+2)//10
    TempFile = np.transpose(np.loadtxt('ParticleSize_'+str(FirstDigit)+'_'+str(LastDigit)+os.sep+'EquilibriumPosition',skiprows=1))
    EquilibriumPoint[iterator]=TempFile[1,-1]

plt.figure(1)

plt.plot(ParticleSize,EquilibriumPoint)

plt.show()