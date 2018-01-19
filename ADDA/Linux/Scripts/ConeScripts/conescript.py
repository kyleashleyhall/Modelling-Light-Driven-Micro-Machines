# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 19:02:12 2017

@author: Marius
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import timeit as tt
import random
import os

def cone(R, h, dipsep):
    V = np.pi*(R**2)*h/3
    dVp = int(V/(dipsep**3))
    if dVp % 2 == 0:
        dV = dVp
    else:
        dV = dVp+1
#    theta = 2*np.arctan(R/h) #not used in calculation.
    dipoles = np.zeros([dV*2,3])
    z = 0
    i = 0
    while z <= h:
        r = R*(h-z)/h
        y = r
        yp = -r
        while yp <= y:
            x = np.sqrt(r**2-yp**2)
            xp = -x
            while xp <= x:
                dipoles[i,0] = xp
                dipoles[i,1] = yp
                dipoles[i,2] = z
                xp += dipsep
                i += 1
                
            yp += dipsep            
        z += dipsep            
    c = 0           
    for e in range(len(dipoles)): #calculate the number of duplicate 0's at the end of the array
        if dipoles[e,0] == 0:
            if dipoles[e,1] == 0:
                if dipoles[e,2] == 0:
                    c += 1

    dipoles = np.round(np.delete(dipoles,np.s_[2*dV-c::],0)/dipsep) #delete all duplicates of 0, divide by dip sep to get integers, then round to nearest number
    return(dipoles)
    
rt = (cone(3,1,0.01))
rt2 = (np.transpose(rt)) #for plotting purposes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(rt2[0], rt2[1], rt2[2])
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
plt.plot(rt2[0],rt2[2], 'g')
plt.show()
print(len(rt))
#c = 0
#for e in range(len(rt)):
#    if rt[e,0]+rt[e,1]+rt[e,2] == 0:
#        c += 1
#print(c)

def prop(l, w, dipsep):
    V = l*(w**2)
    dVp = int(V/(dipsep**3))
    if dVp % 2 == 0:
        dV = dVp
    else:
        dV = dVp+1
    dipoles = np.zeros([50000,3])
    z = 0
    i = 0
    while z <= l:
        r = w*np.sin(np.pi*z/l)
        yp = -r
        while yp <= r:
            xp = -dipsep
            while xp <= dipsep:
                #dipoles[i,0] = xp
                #dipoles[i,1] = yp
                #dipoles[i,2] = z
                xp += dipsep
                i += 1
            yp + dipsep
        z += dipsep
    print(i)
    c = 0
    for e in range(len(dipoles)): #calculate the number of duplicate 0's at the end of the array
        if dipoles[e,0] == 0:
            if dipoles[e,1] == 0:
                if dipoles[e,2] == 0:
                    c += 1

    dipoles = np.round(np.delete(dipoles,np.s_[50000-c::],0)/dipsep) #delete all duplicates of 0, divide by dip sep to get integers, then round to nearest number
    return(dipoles)

#rt = prop(1,3,0.1)
#rt2 = np.transpose([rt])
#plt.plot(rt2[0], rt2[2])
#plt.show()





        
        
    