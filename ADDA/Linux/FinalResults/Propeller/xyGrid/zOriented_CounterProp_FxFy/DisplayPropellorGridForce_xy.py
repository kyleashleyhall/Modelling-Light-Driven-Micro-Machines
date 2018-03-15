# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d
    
TitleString=r"Forces on a z oriented propellor in the xy plane"    


ax = plt.gca()
ax.set_aspect('equal', 'box')
Quiv = np.transpose(np.loadtxt('TotalForce'))
BeamCircle=plt.Circle((0, 0), 1, color='b', fill=False,label='Beam Focus')
ParticleInBeamCircle=plt.Circle((0, 0), 2, color='r', fill=False,label='Particle in Beam')
ForcePlot=plt.quiver(Quiv[0], Quiv[1], Quiv[3], Quiv [4],scale=1e-13,scale_units='x',label=r"Force ($100 fN$)")
plt.title(TitleString)
plt.xlabel(r"Particle position, x ($\mu m$)")
plt.ylabel(r"Particle position, y ($\mu m$)")
ax.add_artist(BeamCircle)
ax.add_artist(ParticleInBeamCircle)
plt.legend(handles=[ForcePlot,BeamCircle,ParticleInBeamCircle])
plt.show()