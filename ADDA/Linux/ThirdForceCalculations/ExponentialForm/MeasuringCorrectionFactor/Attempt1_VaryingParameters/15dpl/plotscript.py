# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 11:36:23 2018

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt

file1 = np.loadtxt('CorrectionFactors', skiprows=1)

#==============================================================================
# plt.plot(file1[:,1],file1[:,-1])
# plt.show()
# plt.plot(file1[:,2],file1[:,-1])
# plt.show()
# plt.plot(file1[:,3],file1[:,-1])
# plt.show()
#==============================================================================
plt.hist(file1[:,-1],100)
plt.show()

for i in range(11):
    t = np.zeros([1])
    for j in range(len(file1)):
        if file1[i,3] == file1[j,3]:
            t = np.append(t,np.array([file1[j,4]]),axis=0)
    t = np.delete(t,[0],axis=0)
    print('Particle size is '+str(file1[i,3]))
    plt.plot(t)
    plt.show()

