# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:30:26 2018

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt

File = np.loadtxt('Positions')

plt.figure(1)
plt.plot(File[:,0], File[:,1])
plt.axhline(y=2),plt.axhline(y=1),plt.axhline(y=-1),plt.axhline(y=-2)
plt.ylim(-2.5, 2.5)
plt.title('Position in y')


plt.figure(2)
plt.hist(File[:,2], bins=20)
plt.title('Brownian Force in y')

plt.show()



