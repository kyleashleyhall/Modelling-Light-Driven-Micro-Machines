# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:34:48 2018

@author: Marius
"""

import numpy as np
import matplotlib.pyplot  as plt
from mpl_toolkits.mplot3d import Axes3D

def Sphere(R, dipsep):
    V = (4/3)*np.pi*(R**3)
    dV = int(V/(dipsep**3))
    dipoles = np.zeros([0,3])
    biggie = int(np.round(R/dipsep)) 
    z,y,x = -biggie,-biggie,-biggie
    while z <= biggie:
        while y <= biggie:
            while x <= biggie:
                if np.sqrt(z**2 + y**2 + x**2) <= biggie:
                    print
                    dipoles = np.append(dipoles,np.array([[x,y,z]]), axis=0)
                x += 1
            x = -biggie
            y += 1
        y = -biggie
        z += 1
    print(len(dipoles))
    return dipoles
#Radius 1, dipsep 0.05 gives 33401 dipoles
#Radius 1 dipsep 0.1 gives 4169 dipoles
DipoleArray = Sphere(1, 0.1)
np.savetxt('spherefile', DipoleArray, fmt='%d')
xyzrows = np.transpose([DipoleArray])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xyzrows[0], xyzrows[1], xyzrows[2])
ax.set_xlabel('0 Column')
ax.set_ylabel('1 Column')
ax.set_zlabel('2 Column')
plt.show()        