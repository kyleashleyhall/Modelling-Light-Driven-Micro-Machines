# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 10:00:30 2018

@author: Marius
"""

import numpy as np
import glob
import os
import subprocess
import time
import shutil
import random
#import matplotlib.pyplot as plt
import scipy.constants as constants

def DragCoef(nu, r):
    return 6*np.pi*nu*r
    
def BrownianForce(Dragcoefficient, tempertature,sigma):
    Boltzmann=1.38064852e-23
    return np.sqrt(2*Dragcoefficient*(Boltzmann)*(tempertature+273))*random.gauss(0,sigma)
    
def DiffusionCoefficient(temperature, viscosity, radius):
    Boltzmann=1.38064852e-23
    return (Boltzmann*(temperature+273)) / (6*np.pi*viscosity*radius)
    
def PositionChange(Force, Dragcoefficient, timestep):
    return ((Force*timestep)/Dragcoefficient)
        
def DipSep(Singleaxis):
    dx = np.zeros([len(Singleaxis)-1])
    for i in range(len(Singleaxis)-1):
        dx[i] = Singleaxis[i+1] - Singleaxis[i]
        if (dx[i]==0):
            dx[i]=np.inf
    return min(np.absolute(dx))

Force = np.transpose(np.loadtxt('TotalForce'))

nu = 8.891e-4
r = 1e-6
Drag_Coefficient = DragCoef(nu,r)
Temperature = 20
D = DiffusionCoefficient(Temperature, nu, r)
arbvalue = 0
sigma = 100
sigma_step = 1
while arbvalue == 0: 
    PPositionArray = np.zeros([0,5])
    t_0 = 0
    t_end = 1
    t_step = 1e-4
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
    if -1e-15 < (gradient - (2*D)) < 1e-15:
        arbvalue += 1
        
    elif sigma_step <= 1e-10:
        arbvalue += 1
        
    elif gradient < (2*D):
        sigma += sigma_step
        
    elif gradient > (2*D):
        sigma -= sigma_step
        sigma_step *= 0.5
        sigma += sigma_step
        
print(sigma)

#Preliminary Dynamic variables
t_0 = 0
t_end = 50
t_step = 1e-4
randomval = np.random.random()

#Particle start point start point
x_position, y_position, z_position = 0,0,0

PositionArray = np.array([[t_0,y_position,0]])

while t_0 <= t_end:
    Force_y = Force[1,np.argmin(np.abs(y_position - Force[0]))]
    Brownian = BrownianForce(Drag_Coefficient, Temperature, sigma)
    Force_Total = Force_y + Brownian    
    y_position += PositionChange(Force_Total, Drag_Coefficient, t_step)*(1e6)
    t_0 += t_step    
    PositionArray = np.append(PositionArray, np.array([[t_0, y_position, Brownian]]), axis=0)
    np.savetxt(('Positions'+str(randomval)), PositionArray, fmt='%e', delimiter='  ')