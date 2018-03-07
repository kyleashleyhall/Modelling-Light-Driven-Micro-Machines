# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d

BeamWidth=int(((os.getcwd())[-4:-3]))+((int(((os.getcwd())[-2:]))/100))
TitleString=r"Torque about the z axis on a cone in the xy plane for a "+str(BeamWidth)+" $\mu m$ beam width"

ax = plt.gca()
ax.set_aspect('equal', 'box')
TorqueGrid = np.loadtxt('zTorqueGrid')
xValues=np.loadtxt('xPositions')
yValues=np.loadtxt('yPositions')
BeamCircle=plt.Circle((0, 0), BeamWidth, color='b', fill=False,label='Beam Focus')
ParticleInBeamCircle=plt.Circle((0, 0), (BeamWidth+1), color='r', fill=False,label='Particle in Beam')
TorqueGrid*=1e-6
TorquePlot=plt.pcolormesh(xValues,yValues,TorqueGrid,cmap='bwr',vmin=-1e-20,vmax=1e-20)
plt.title(TitleString)
plt.xlabel(r"Particle position, x ($\mu m$)")
plt.ylabel(r"Particle position, y ($\mu m$)")
ax.add_artist(BeamCircle)
ax.add_artist(ParticleInBeamCircle)
plt.legend(handles=[TorquePlot,BeamCircle,ParticleInBeamCircle])
Colorbar=plt.colorbar()
Colorbar.set_label("Torque about z ($Nm$)")
plt.show()