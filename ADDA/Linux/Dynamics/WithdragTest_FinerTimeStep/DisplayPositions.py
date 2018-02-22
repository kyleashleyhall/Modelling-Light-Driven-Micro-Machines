# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Positions = np.transpose(np.loadtxt('ParticlePositions'))

plt.figure(1)

plt.subplot(311)
plt.plot(Positions[0],Positions[1])

plt.subplot(312)
plt.plot(Positions[0],Positions[2])

plt.subplot(313)
plt.plot(Positions[0],Positions[3])

plt.show()