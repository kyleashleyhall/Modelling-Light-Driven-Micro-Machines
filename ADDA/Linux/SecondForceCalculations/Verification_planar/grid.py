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

def DipSep(Singleaxis):
    dx = np.zeros([len(Singleaxis)-1])
    for i in range(len(Singleaxis)-1):
        dx[i] = Singleaxis[i+1] - Singleaxis[i]
        if (dx[i]==0):
            dx[i]=np.inf
    return min(np.absolute(dx))

def FileSlice(Fname):
	filename = np.loadtxt(Fname, skiprows=1)	
	x,y,z = filename[:,0], filename[:,1], filename[:,2]
	varX = filename[:,4] + 1j*filename[:,5]
	varY = filename[:,6] + 1j*filename[:,7]
	varZ = filename[:,8] + 1j*filename[:,9]
	return np.vstack([x,y,z]), np.vstack([varX,varY,varZ])

def Grid(PositionsArray,GridType):
    Sep = DipSep(PositionsArray[0,:]) #dipsep only works in x axis
    GridSize = np.zeros([3], dtype=int) #Dimensions of Grid in x,y,z
    for i in range(3):
        minLocation = min(IntFieldRawPositions[i,:])
        maxLocation = max(IntFieldRawPositions[i,:])
        GridSize[i] = np.round((maxLocation-minLocation)/Sep)+1
    if (GridType=='float'):
        grid = np.zeros([GridSize[0],GridSize[1],GridSize[2],3], dtype=float)
        grid[:] = np.nan
    elif (GridType=='complex'):
        grid = np.zeros([GridSize[0],GridSize[1],GridSize[2],3], dtype=complex)
        grid[:] = np.nan
        
    return grid


IntFieldRawPositions = FileSlice('IntField-Y')[0]
IntFieldRawValues = FileSlice('IntField-Y')[1]
IntFieldGridPositions = Grid(IntFieldRawPositions,'float')
#print(IntFieldGridPositions)
IntFieldGridValues = Grid(IntFieldRawPositions,'complex')

min_xLocation=min(IntFieldRawPositions[0,:]) #Min x
min_yLocation=min(IntFieldRawPositions[1,:]) #Min y
min_zLocation=min(IntFieldRawPositions[2,:]) #Min z
DipoleSep=DipSep(IntFieldRawPositions[0,:]) #The Dipole seperation

for i in range(len(IntFieldRawPositions[1])):
    xLocation=np.int(np.round((IntFieldRawPositions[0,i]-min_xLocation)/DipoleSep)) #Find the x Location in the grid
    yLocation=np.int(np.round((IntFieldRawPositions[1,i]-min_yLocation)/DipoleSep)) #Find the y Location in the grid
    zLocation=np.int(np.round((IntFieldRawPositions[2,i]-min_zLocation)/DipoleSep)) #Find the z Location in the grid
    IntFieldGridPositions[xLocation,yLocation,zLocation,0]=IntFieldRawPositions[0,i]
    IntFieldGridPositions[xLocation,yLocation,zLocation,1]=IntFieldRawPositions[1,i]
    IntFieldGridPositions[xLocation,yLocation,zLocation,2]=IntFieldRawPositions[2,i]
    IntFieldGridValues[xLocation,yLocation,zLocation,0]=IntFieldRawValues[0,i]
    IntFieldGridValues[xLocation,yLocation,zLocation,1]=IntFieldRawValues[1,i]
    IntFieldGridValues[xLocation,yLocation,zLocation,2]=IntFieldRawValues[2,i]