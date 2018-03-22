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
import matplotlib.pyplot as plt

def DragCoef(nu, r):
    return 6*np.pi*nu*r
    
def DiffusionCoefficient(temperature, viscosity, radius):
    Boltzmann=1.38064852e-23
    return (Boltzmann*(temperature+273)) / (6*np.pi*viscosity*radius)
    
def BrownianForce(Dragcoefficient,temperature, sigma):
    Boltzmann=1.38064852e-23
    return np.sqrt(2*Dragcoefficient*(Boltzmann)*(temperature+273))*random.gauss(0,sigma)

    
def PositionChange(Force, Dragcoefficient, timestep):
    return ((Force*timestep)/Dragcoefficient)
        
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
r = 3e-7

#Preliminary Dynamic variables
t_end = 1
t_step = 1e-4





#Array to track particle position
StartTime=time.clock()

Drag_Coefficient = DragCoef(nu,r)
D = DiffusionCoefficient(Temperature, nu, r)
arbvalue = 0
sigma = 100
sigma_step = 1
while arbvalue == 0: 
    PPositionArray = np.zeros([0,5])
    t_0 = 0
    #Particle start point start point
    x_position = 0
    x_positionsquare = 0
    
    while t_0 < t_end:      
        
        #Generate the Brownian "Force"
        B_x = BrownianForce(Drag_Coefficient, Temperature, sigma)

        EstimatedParticleForce3 = np.array([[B_x]])
                    
        #Calculate the change in position of the particle
        x = PositionChange(EstimatedParticleForce3[0], Drag_Coefficient, t_step)
        x_position += x[0]
        x_positionsquare += x[0]**2
        t_0 += t_step  
        PPositionArray = np.append(PPositionArray, np.array([[t_0, D, x_position, x_positionsquare, EstimatedParticleForce3[0]]]),axis=0)
    
    gradient = np.mean(np.gradient(PPositionArray[:,3], t_step))
    if gradient < (2*D):
        sigma += sigma_step
        
    if gradient > (2*D):
        sigma -= sigma_step
        sigma_step *= 0.5
        sigma += sigma_step

    if -1e-15 < (gradient - (2*D)) < 1e-15:
        arbvalue += 1
        
    if sigma_step <= 1e-10:
        arbvalue += 1
        print('This value of sigma is accurate to 9 decimal places')
        
    EndTime=time.clock()
    TimeRecordings=np.array([[(EndTime-StartTime)]])        
    np.savetxt('ParticlePositions', PPositionArray, fmt='%e', delimiter=' ')
    TimeLogPath = str(os.getcwd())+str(os.sep+'TimeLog')	
    with open(TimeLogPath, 'wb') as f:
        f.write(b'Time(s)\n')
        np.savetxt(f, TimeRecordings, fmt='%.10f', delimiter=' ')
        
    File = np.loadtxt('ParticlePositions', skiprows=1)
    plt.plot(File[:,0], File[:,3], 'k', label='mean square displacement')
    plt.plot(File[:,0], (2*File[:,0]*File[:,1]), 'b', label='2Dt')
    plt.legend(loc='upper left')
    plt.show()
    plt.plot(File[:,0], File[:,2])
    plt.show()
        
        
    



        


