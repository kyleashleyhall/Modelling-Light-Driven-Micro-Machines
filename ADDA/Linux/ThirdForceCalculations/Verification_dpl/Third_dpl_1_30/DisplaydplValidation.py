# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

FileName=input('Filename: ')

Forces = np.transpose(np.loadtxt(FileName, skiprows=1))

plt.figure(1)

plt.suptitle('dpl Validation Results', fontsize=16)

plt.subplot(321)
plt.plot(Forces[0],Forces[1])
plt.xlabel('Dipoles per Lambda')
plt.xlabel('Dipoles per Lambda')

plt.subplot(322)
plt.title("Column 4")
plt.plot(Forces[0],Forces[4])
plt.xlabel('Dipoles per Lambda')

plt.subplot(323)
plt.plot(Forces[0],Forces[2])
plt.xlabel('Dipoles per Lambda')

plt.subplot(324)
plt.title("Column 5")
plt.plot(Forces[0],Forces[5])
plt.xlabel('Dipoles per Lambda')

plt.subplot(325)
plt.plot(Forces[0],Forces[3])
plt.xlabel('Dipoles per Lambda')

plt.subplot(326)
plt.title("Column 6")
plt.plot(Forces[0],Forces[6])
plt.xlabel('Dipoles per Lambda')

plt.show()