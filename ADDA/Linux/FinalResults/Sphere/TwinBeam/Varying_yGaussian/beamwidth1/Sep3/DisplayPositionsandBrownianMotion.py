# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:30:26 2018

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt

File = np.transpose(np.loadtxt('Positions'))

plt.figure(1)
plt.plot(File[:,0],File[:,1])

plt.show()

plt.figure(2)
plt.plot(File[:,0],File[:,2])

plt.show()



