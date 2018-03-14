# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 20:48:12 2018

@author: Marius
"""

import os
import numpy as np
import glob
import matplotlib.pyplot as plt

File = np.loadtxt('EquilibriumPosition', skiprows=1)

plt.plot(File[:,1], File[:,2])
plt.xlabel('Position in z')
plt.ylabel('Force in z')
plt.show()