# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:30:26 2018

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt

File = np.loadtxt('ParticlePositions')

plt.figure(1)
plt.plot(File[:,0], File[:,1])
plt.title('Position in y')

plt.figure(2)
plt.title('Brownian Force in y')
plt.hist(File[:,3], bins=20)

plt.figure(3)
plt.hist(File[:,2])
plt.title('Force in y')

plt.show()



