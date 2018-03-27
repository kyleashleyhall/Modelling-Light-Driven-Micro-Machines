# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 15:03:01 2018

@author: Marius and Kyle
"""

import numpy as np
import os
import time
import random
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
    
def ADDAForceConversion(Force,ElectricFieldStrength):

	Force=Force*(((ElectricFieldStrength)**2)/((constants.c)**2))*(10**(-5))
	
	return(Force)
	
def OurForceConversion(Force,CorrectionFactor,ElectricFieldStrength):
    
   return(ADDAForceConversion((Force*CorrectionFactor),ElectricFieldStrength))
   
def ElectricFieldStrengthCalc(DielConstant,Power,BeamWidth): #Power in Watts, BeamWidth in micro m
    
    Impedence=((constants.mu_0)/((constants.epsilon_0)*(DielConstant)))**0.5
    ElectricFieldStrength=((Impedence*Power)/((constants.pi)*(((BeamWidth*(1e-6))/2)**2)))**0.5
    
    return(ElectricFieldStrength)

#Preliminary Variables	
CD = os.getcwd()
size = 0.6
BeamWidth = 1 #In micro m
Temperature = 20 #Degrees C
nu = 8.891e-4
r = (size/2e6)

Drag_Coefficient = DragCoef(nu,r)
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


Power = 5e-3 #In Watts
MediumDielectricConstant=87.740-(0.40008*Temperature)+(9.398e-4*(Temperature**2))-(1.410e-6*(Temperature**3))
ElectricFieldStrength=ElectricFieldStrengthCalc(MediumDielectricConstant,Power,BeamWidth) #V/m


#Preliminary Dynamic variables
t_0 = 0
t_end = 60
t_step = 1e-4

#Particle start point start point
x_position, y_position, z_position = 0,0,0


#Load Force Field

ForcesFile=np.transpose(np.loadtxt('ForcesY'))
Positions=ForcesFile[0]
Forces=ForcesFile[1]


#Array to track particle position
PPositionArray = np.zeros([1,4])
StartTime=time.clock()
while t_0 <= t_end:


    #Generate the Brownian "Force"
    Drag_Coefficient = DragCoef(nu,r)
    B_y = BrownianForce(Drag_Coefficient, Temperature, sigma)
    
    #Calculate Estimated Particle Force
    if np.any((y_position<Positions)):
        if(np.all((y_position<Positions))):
            EstimatedParticleForce=0 #Outside Grid
            
        else:
            PositionArgumentAbove=np.argwhere(y_position<Positions)[0,0]
            ForceAbove=Forces[PositionArgumentAbove]
            PositionArgumentBelow=PositionArgumentAbove-1
            ForceBelow=Forces[PositionArgumentBelow]
            ForceDifference=ForceAbove-ForceBelow
            Seperation=Positions[PositionArgumentAbove]-Positions[PositionArgumentBelow]
            EstimatedParticleForce=ForceBelow+(((y_position-Positions[PositionArgumentBelow])/Seperation)*ForceDifference)
            
            
    else:
        EstimatedParticleForce=0 #Outside Grid
    
    EstimatedParticleForce = EstimatedParticleForce+B_y
                    
    #Calculate the change in position of the particle
    y = PositionChange(EstimatedParticleForce, Drag_Coefficient, t_step)*(1e6)
    y_position += y
    t_0 += t_step
    PPositionArray = np.append(PPositionArray, np.array([[t_0,y_position,EstimatedParticleForce, B_y]]),axis=0)
    EndTime=time.clock()
    
np.savetxt('ParticlePositions', PPositionArray, fmt='%e', delimiter=' ')
print(EndTime-StartTime)