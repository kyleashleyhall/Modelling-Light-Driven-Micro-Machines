# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Forces = np.transpose(np.loadtxt('ParticlePositions'))

plt.figure(1)

plt.subplot(311)
plt.title('Particle Position')
plt.plot(Forces[0],Forces[3])

plt.subplot(312)
plt.title('Particle Forces')
plt.plot(Forces[0],Forces[6])

plt.subplot(313)
plt.title('Brownian Force')
plt.hist(Forces[7])

plt.show()