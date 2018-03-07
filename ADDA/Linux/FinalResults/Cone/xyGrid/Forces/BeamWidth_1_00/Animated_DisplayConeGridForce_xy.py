# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import animation

BeamWidth=int(((os.getcwd())[-4:-3]))+((int(((os.getcwd())[-2:]))/100))
TitleString=r"Forces on a cone in the xy plane for a "+str(BeamWidth)+" $\mu m$ beam width"

fig,ax = plt.subplots(1,1)
ax.set_aspect('equal', 'box')

Quiv_Read = np.transpose(np.loadtxt('BeamWidth_0_20'+os.sep()+'TotalForce'))
Quiv=np.zeros([[17,np.size(Quiv_Read,axis=0),np.size(Quiv_Read,axis=1)]])
BeamWidths=np.zeros(17)
Quiv[0,:,:]=Quiv_Read
BeamWidths[0]=0.2

Quiv_Read = np.transpose(np.loadtxt('BeamWidth_0_5'+os.sep()+'TotalForce'))

BeamCircle=plt.Circle((0, 0), BeamWidth, color='b', fill=False,label='Beam Focus')
ParticleInBeamCircle=plt.Circle((0, 0), (BeamWidth+1), color='r', fill=False,label='Particle in Beam')
ForcePlot=plt.quiver(Quiv[0,0], Quiv[0,1], Quiv[0,3], Quiv[0,4],scale=1e-13,scale_units='x',label=r"Force ($100 fN$)")
plt.title(TitleString)
plt.xlabel(r"Particle position, x ($\mu m$)")
plt.ylabel(r"Particle position, y ($\mu m$)")
ax.add_artist(BeamCircle)
ax.add_artist(ParticleInBeamCircle)
plt.legend(handles=[ForcePlot,BeamCircle,ParticleInBeamCircle])

def update_quiver(num,ForcePlot):
    
    print(num)
    
    return(ForcePlot)
    
anim = animation.FuncAnimation(fig, update_quiver, fargs=(ForcePlot),interval=50, blit=False)

plt.show()