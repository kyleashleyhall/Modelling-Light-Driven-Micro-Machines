# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d
    
TitleString=r"Torque about the x axis on a y oriented propellor in the xy plane"

ax = plt.gca()
ax.set_aspect('equal', 'box')
TorqueGrid = np.loadtxt('xTorqueGrid')
xValues=np.loadtxt('xPositions')
yValues=np.loadtxt('yPositions')
BeamCircle=plt.Circle((0, 0), 1, color='b', fill=False,label='Beam Focus')
ParticleInBeamCircle=plt.Circle((0, 0), 2, color='r', fill=False,label='Particle in Beam')
TorqueGrid*=1e-6
TorquePlot=plt.pcolormesh(xValues,yValues,TorqueGrid,cmap='bwr')
plt.title(TitleString)
plt.xlabel(r"Particle position, x ($\mu m$)")
plt.ylabel(r"Particle position, y ($\mu m$)")
ax.add_artist(BeamCircle)
ax.add_artist(ParticleInBeamCircle)
plt.legend(handles=[TorquePlot,BeamCircle,ParticleInBeamCircle])
Colorbar=plt.colorbar()
Colorbar.set_label("Torque about x ($Nm$)")
plt.show()