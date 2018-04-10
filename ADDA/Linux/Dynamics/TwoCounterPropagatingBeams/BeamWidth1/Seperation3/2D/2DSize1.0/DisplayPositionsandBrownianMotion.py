# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:30:26 2018

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt

File = np.loadtxt('ParticlePositions')

plt.figure(1)
plt.subplot(311)
plt.plot(File[:,0], File[:,1])
plt.title('Position in x')

plt.subplot(312)
plt.plot(File[:,0], File[:,2])
plt.title('Position in y')

plt.subplot(313)
plt.plot(File[:,0], np.sqrt(File[:,1]**2+File[:,2]**2))
plt.title('Radial distance from the centre')

plt.figure(2)
plt.subplot(311)
plt.title('Brownian Force in x')
plt.hist(File[:,7], bins=20)

plt.subplot(312)
plt.title('Brownian Force in y')
plt.hist(File[:,8], bins=20)

plt.subplot(313)
plt.title('Brownian Force in z')
plt.hist(File[:,9], bins=20)

plt.figure(3)
plt.subplot(311)
plt.hist(File[:,4])
plt.title('Force in x')

plt.subplot(312)
plt.hist(File[:,5])
plt.title('Force in y')

plt.subplot(313)
plt.hist(File[:,6])
plt.title('Force in z')

plt.show()



