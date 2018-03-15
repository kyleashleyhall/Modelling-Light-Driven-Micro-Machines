# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 12:45:10 2018

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt

File = np.loadtxt('ParticlePositions', skiprows=1)
plt.plot(File[:,0], File[:,3], 'k', label='mean square displacement')
plt.plot(File[:,0], (2*File[:,0]*File[:,1]), 'b', label='2Dt')
plt.legend(loc='upper left')
plt.show()