# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d

Quiv = np.transpose(np.loadtxt('CalculatedForces',skiprows=1))

print(sum(Quiv[4]))
print(sum(Quiv[5]))
print(sum(Quiv[6]))

plt.quiver(Quiv[1], Quiv[2], Quiv[5], Quiv [6])
plt.show()