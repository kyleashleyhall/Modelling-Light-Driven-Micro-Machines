# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Forces1 = np.transpose(np.loadtxt('ForcesY'))
PotentialX = np.gradient(Forces1[1])
PotentialY = np.gradient(Forces1[2])
PotentialZ = np.gradient(Forces1[3])
plt.figure(1)

plt.subplot(311)
plt.plot(Forces1[0],PotentialX)

plt.subplot(312)
plt.plot(Forces1[0],PotentialY)

plt.subplot(313)
plt.plot(Forces1[0],PotentialZ)

plt.show()