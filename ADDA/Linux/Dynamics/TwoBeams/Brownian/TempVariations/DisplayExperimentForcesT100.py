# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Forces = np.transpose(np.loadtxt('ParticlePositionsT100'))

plt.figure(1)
plt.title('Positions in time')
plt.subplot(311)
plt.plot(Forces[0],Forces[1])

plt.subplot(312)
plt.plot(Forces[0],Forces[2])

plt.subplot(313)
plt.plot(Forces[0],Forces[3])

plt.show()

plt.figure(1)
plt.title('Forces in Time')
plt.subplot(311)
plt.plot(Forces[0],Forces[4])

plt.subplot(312)
plt.plot(Forces[0],Forces[5])

plt.subplot(313)
plt.plot(Forces[0],Forces[6])

plt.show()