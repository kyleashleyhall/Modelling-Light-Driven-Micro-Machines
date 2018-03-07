# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Forces = np.transpose(np.loadtxt('ParticlePositions'))


plt.figure(1)
plt.title('Forces in Time')
plt.subplot(211)
plt.plot(Forces[0],Forces[3])

plt.subplot(212)
plt.plot(Forces[0],Forces[4])


plt.figure(2)
plt.title('Positions in time')
plt.subplot(211)
plt.plot(Forces[0],Forces[1])

plt.subplot(212)
plt.plot(Forces[0],Forces[2])

plt.show()
