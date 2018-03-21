# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d
    
TitleString=r"Force in the y axis for a sphere in a varying twin beam arrangement"

ax = plt.gca()
ax.set_aspect('equal', 'box')
TorqueGrid = np.loadtxt('yForceGrid')
xValues=np.loadtxt('xAxis')
yValues=np.loadtxt('yAxis')
TorquePlot=plt.pcolormesh(xValues,yValues,TorqueGrid,cmap='bwr')
plt.title(TitleString)
plt.xlabel(r"Particle position, y ($\mu m$)")
plt.ylabel(r"Beam seperation ($\mu m$)")
Colorbar=plt.colorbar()
Colorbar.set_label("Force in y ($N$)")
plt.show()