# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

def shitIntegral(val1, val2, increment):
    return ((val2+val1)/2)*increment

Forces = np.transpose(np.loadtxt('TotalForce'))

#==============================================================================
# plt.figure(1)
# 
# plt.subplot(311)
# plt.plot(Forces[1],Forces[3])
# 
# plt.subplot(312)
# plt.plot(Forces[1],Forces[4])
# 
# plt.subplot(313)
# plt.plot(Forces[1],Forces[5])
# 
# plt.show()
#==============================================================================
#Potential = np.zeros(len(Forces[1])-1)
increment = Forces[1,1] - Forces[1,0]
Potential = np.gradient(Forces[4])
#==============================================================================
# for i in range(len(Potential)):
#     Potential[i] = -np.nan_to_num(1/shitIntegral(Forces[5,i], Forces[5,i+1], increment))
# print(Potential)
#==============================================================================
y = np.linspace(min(Forces[1]), max(Forces[1]), (len(Forces[1])-1))
plt.plot(Forces[1], Potential)
plt.show()