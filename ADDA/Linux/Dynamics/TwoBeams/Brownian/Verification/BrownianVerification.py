# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 15:03:01 2018

@author: Marius and Kyle
"""

import numpy as np
import glob
import os
import subprocess
import time
import shutil
import random
import scipy.constants as constants

def DragCoef(nu, r):
    return 6*np.pi*nu*r
    
def DiffusionCoefficient(temperature, viscosity, radius):
    Boltzmann=1.38064852e-23
    return (Boltzmann*(temperature+273)) / (6*np.pi*viscosity*radius)
    
def BrownianForce(DiffusionCoef):
    return np.sqrt(2*DiffusionCoef)*random.gauss(0,1)
    
def PositionChange(Force, Dragcoefficient, timestep):
    return ((Force*timestep)/Dragcoefficient)*(1e-6)
        
def DipSep(Singleaxis):
    dx = np.zeros([len(Singleaxis)-1])
    for i in range(len(Singleaxis)-1):
        dx[i] = Singleaxis[i+1] - Singleaxis[i]
        if (dx[i]==0):
            dx[i]=np.inf
    return min(np.absolute(dx))
    


#Preliminary Variables	

BeamWidth=1 #In micro m
Temperature=20 #Degrees C


nu = 8.891e-4
r = 1e-6

#Preliminary Dynamic variables
t_0 = 0
t_end = 1
t_step = 1e-4



#Particle start point start point
x_position = 0
x_positionsquare = 0

#Array to track particle position
PPositionArray = np.zeros([1,5])
StartTime=time.clock()
while t_0 < t_end:
    #Generate the Brownian "Force"
    Drag_Coefficient = DragCoef(nu,r)
    D = DiffusionCoefficient(Temperature, nu, r)
    B_x = BrownianForce(D)

    EstimatedParticleForce3 = np.array([[B_x]])
                    
    #Calculate the change in position of the particle
    x = PositionChange(EstimatedParticleForce3[0], Drag_Coefficient, t_step)
    x_position += x[0]
    x_positionsquare += x[0]**2
    t_0 += t_step  
    PPositionArray = np.append(PPositionArray, np.array([[t_0, D, x_position, x_positionsquare, EstimatedParticleForce3[0]]]),axis=0)
         



        
EndTime=time.clock()
TimeRecordings=np.array([[(EndTime-StartTime)]])        
np.savetxt('ParticlePositions', PPositionArray, fmt='%e', delimiter=' ')
TimeLogPath = str(os.getcwd())+str(os.sep+'TimeLog')	
with open(TimeLogPath, 'wb') as f:
    f.write(b'Time(s)\n')
    np.savetxt(f, TimeRecordings, fmt='%.10f', delimiter=' ')

