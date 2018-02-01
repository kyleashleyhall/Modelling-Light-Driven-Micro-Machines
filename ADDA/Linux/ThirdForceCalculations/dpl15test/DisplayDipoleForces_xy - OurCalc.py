# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d

Quiv = np.transpose(np.loadtxt('CalculatedForces'))

print(sum(Quiv[3]))
print(sum(Quiv[4]))
print(sum(Quiv[5]))

plt.quiver(Quiv[0], Quiv[1], Quiv[3], Quiv [4])
plt.show()