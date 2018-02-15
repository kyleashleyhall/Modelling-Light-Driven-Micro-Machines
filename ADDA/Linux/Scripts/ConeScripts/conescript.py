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

    dipoles = np.int16(np.round(np.delete(dipoles,np.s_[2*dV-c::],0)/dipsep)) #delete all duplicates of 0, divide by dip sep to get integers, then round to nearest number
    print(len(dipoles))
    return(dipoles)
    

#np.savetxt('conefile', rt, fmt='%d')


def prop(l, w, dipsep):
    V = l*(w**2)
    dVp = int(V/(dipsep**3))
    if dVp % 2 == 0:
        dV = dVp
    else:
        dV = dVp+1
    dipoles1 = np.zeros([dV,3])
    z = 0
    i = 0
    while z <= l:
        r = w
        yp = 0
        while yp <= r:
            xp = yp
            while xp <= r:
                dipoles1[i,0] = z
                dipoles1[i,1] = yp
                dipoles1[i,2] = xp
                xp += dipsep
                i += 1
            yp += dipsep
        z += dipsep
    c = 0
    for e in range(len(dipoles1)): #calculate the number of duplicate 0's at the end of the array
        if dipoles1[e,0] == 0:
            if dipoles1[e,1] == 0:
                if dipoles1[e,2] == 0:
                    c += 1

    dipoles1 = np.round(np.delete(dipoles1,np.s_[dV-c::],0)/dipsep)
    #delete all duplicates of 0, divide by dip sep to get integers, then round to nearest number
    dipoles2 = np.zeros([dV,3])
    z = 0
    j = 0
    while z <= l:
        r = w
        yp = 0
        while yp <= r:
            xp = yp
            while xp <= r:
                dipoles2[j,0] = -z
                dipoles2[j,1] = -yp
                dipoles2[j,2] = xp
                xp += dipsep
                j += 1
            yp += dipsep
        z += dipsep
    c = 0
    for e in range(len(dipoles2)): #calculate the number of duplicate 0's at the end of the array
        if dipoles2[e,0] == 0:
            if dipoles2[e,1] == 0:
                if dipoles2[e,2] == 0:
                    c += 1

    dipoles2 = np.round(np.delete(dipoles2,np.s_[dV-c::],0)/dipsep)
    
    dipoles3 = np.zeros([dV,3])
    z = 0
    k = 0
    while z <= l:
        r = w
        yp = 0
        while yp <= r:
            xp = yp
            while xp <= r:
                dipoles3[k,0] = yp
                dipoles3[k,1] = -z
                dipoles3[k,2] = xp
                xp += dipsep
                k += 1
            yp += dipsep
        z += dipsep

    c = 0
    for e in range(len(dipoles3)): #calculate the number of duplicate 0's at the end of the array
        if dipoles3[e,0] == 0:
            if dipoles3[e,1] == 0:
                if dipoles3[e,2] == 0:
                    c += 1

    dipoles3 = np.round(np.delete(dipoles3,np.s_[dV-c::],0)/dipsep)
    
    dipoles4 = np.zeros([dV,3])
    z = 0
    h = 0
    while z <= l:
        r = w
        yp = 0
        while yp <= r:
            xp = yp
            while xp <= r:
                dipoles4[h,0] = -yp
                dipoles4[h,1] = z
                dipoles4[h,2] = xp
                xp += dipsep
                h += 1
            yp += dipsep
        z += dipsep
    c = 0
    for e in range(len(dipoles4)): #calculate the number of duplicate 0's at the end of the array
        if dipoles4[e,0] == 0:
            if dipoles4[e,1] == 0:
                if dipoles4[e,2] == 0:
                    c += 1

    dipoles4 = np.round(np.delete(dipoles4,np.s_[dV-c::],0)/dipsep)
    print(len(dipoles1)+len(dipoles2)+len(dipoles3)+len(dipoles4))
    dipoles = dipoles1
    dipoles = np.append(dipoles,dipoles2,axis=0)
    dipoles = np.append(dipoles,dipoles3,axis=0)
    dipoles = np.append(dipoles,dipoles4,axis=0)
    dipoles = np.vstack({tuple(row) for row in dipoles})
    return dipoles



dipoleSeperation= 0.09955405125-0.03318468375
ConeRadius=1 #1 micro m
ConeHeight=2 #2 micro m

#CONE FILE
dipolearray = cone(ConeRadius,ConeHeight,dipoleSeperation) #cone, with dimensions and increment
np.savetxt('conefile', dipolearray, fmt='%d')
xyzrows = np.transpose([dipolearray])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xyzrows[0], xyzrows[1], xyzrows[2])
ax.set_xlabel('0 Column')
ax.set_ylabel('1 Column')
ax.set_zlabel('2 Column')
plt.show()

#PROPELLOR FILE
dipolearray2 = prop(2,2,dipoleSeperation)
np.savetxt('propellorfile', dipolearray2, fmt='%d')
xyzrows = np.transpose([dipolearray2])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xyzrows[0], xyzrows[1], -xyzrows[2])
ax.set_xlabel('0 Column')
ax.set_ylabel('1 Column')
ax.set_zlabel('2 Column')
plt.show()

        
        
    