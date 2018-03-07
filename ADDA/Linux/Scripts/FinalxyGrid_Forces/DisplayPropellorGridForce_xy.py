# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d

try:
    Alpha=int((os.getcwd())[-11:-9])
except:
    Alpha=int((os.getcwd())[-10:-9])
    
AlphaNumerator=int(Alpha)
AlphaDenominator=40
    
if(Alpha%2==0):
    AlphaNumerator=int(Alpha/2)
    AlphaDenominator=20
    
if(Alpha%4==0):
    AlphaNumerator=int(Alpha/4)
    AlphaDenominator=10
    
if(Alpha%5==0):
    AlphaNumerator=int(Alpha/5)
    AlphaDenominator=8
    
if(Alpha%8==0):
    AlphaNumerator=int(Alpha/8)
    AlphaDenominator=5
    
if(Alpha%10==0):
    AlphaNumerator=int(Alpha/10)
    AlphaDenominator=4
    
if(Alpha%20==0):
    AlphaNumerator=int(Alpha/20)
    AlphaDenominator=2
    
if(AlphaNumerator==1):
    TitleString=r"Forces on a propellor in the xy plane for a rotation of $\frac{\pi}{"+str(AlphaDenominator)+"}$ radians about the z axis"
    
elif(AlphaNumerator==0):
    TitleString=r"Forces on a propellor in the xy plane for a rotation of $0$ radians about the z axis"
    
else:
    TitleString=r"Forces on a propellor in the xy plane for a rotation of $\frac{"+str(AlphaNumerator)+"\pi}{"+str(AlphaDenominator)+"}$ radians about the z axis"    


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