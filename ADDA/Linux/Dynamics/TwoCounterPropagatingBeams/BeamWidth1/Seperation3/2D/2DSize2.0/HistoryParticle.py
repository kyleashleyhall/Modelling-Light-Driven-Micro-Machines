# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:30:26 2018

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt

File = np.loadtxt('ParticlePositions')

plt.figure(1)
plt.plot(File[:,1], File[:,2])
plt.title('Position in x')
plt.ylabel("Y position")


plt.show()



