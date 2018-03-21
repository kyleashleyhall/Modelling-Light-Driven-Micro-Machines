# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

Forces = np.transpose(np.loadtxt('TotalForce'))

plt.figure(1)

plt.plot(Forces[3],Forces[7])

plt.show()