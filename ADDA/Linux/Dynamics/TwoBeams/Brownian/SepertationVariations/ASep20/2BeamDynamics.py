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
    
def BrownianForce(Dragcoefficient, tempertature):
    Boltzmann=1.38064852e-23
    return np.sqrt(2*Dragcoefficient*(Boltzmann)*(tempertature+273))*random.gauss(0,1)
    
def PositionChange(Force, Dragcoefficient, timestep):
    return ((Force*timestep)/Dragcoefficient)*(1e6)
        
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
    
def Forces(DipPol,IntField,IncBeam,DipoleSep):
    
    Force=np.zeros([len(DipPol[0]),7])
    for i in range(len(DipPol[0])):
        Force[i,0]=DipPol[0,i] #Assign the ith x variable from DipPol to the Force file
        Force[i,1]=DipPol[1,i] #Assign the ith y variable from DipPol to the Force file
        Force[i,2]=DipPol[2,i] #Assign the ith z variable from DipPol to the Force file
        
        #Calculate the force in x direction
        linePositionBelow_IncBeam=np.argwhere(((DipPol[0,i]-((3*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]-((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionCurrent_IncBeam=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionAbove_IncBeam=np.argwhere(((DipPol[0,i]+((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((3*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionBelow_IntField=np.argwhere(((DipPol[0,i]-((3*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]-((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionCurrent_IntField=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionAbove_IntField=np.argwhere(((DipPol[0,i]+((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((3*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        CentralDifference=np.ones([2],dtype=bool) #To check if lower and upper part of Central difference is available
        if((len(linePositionBelow_IncBeam)==0)or(len(linePositionBelow_IntField)==0)): #Is there a value below?
            CentralDifference[0]=False
        if((len(linePositionAbove_IncBeam)==0)or(len(linePositionAbove_IntField)==0)): #Is there a value above?
            CentralDifference[1]=False
        if ((CentralDifference[0]==True)&(CentralDifference[1]==True)):
            dEx_dx=((((IncBeam[4,linePositionAbove_IncBeam]+IntField[4,linePositionAbove_IntField])-(IncBeam[4,linePositionBelow_IncBeam]+IntField[4,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[5,linePositionAbove_IncBeam]+IntField[5,linePositionAbove_IntField])-(IncBeam[5,linePositionBelow_IncBeam]+IntField[5,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            dEy_dx=((((IncBeam[6,linePositionAbove_IncBeam]+IntField[6,linePositionAbove_IntField])-(IncBeam[6,linePositionBelow_IncBeam]+IntField[6,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[7,linePositionAbove_IncBeam]+IntField[7,linePositionAbove_IntField])-(IncBeam[7,linePositionBelow_IncBeam]+IntField[7,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            dEz_dx=((((IncBeam[8,linePositionAbove_IncBeam]+IntField[8,linePositionAbove_IntField])-(IncBeam[8,linePositionBelow_IncBeam]+IntField[8,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[9,linePositionAbove_IncBeam]+IntField[9,linePositionAbove_IntField])-(IncBeam[9,linePositionBelow_IncBeam]+IntField[9,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            Force[i,4]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dx)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dx)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dx))
        if ((CentralDifference[0]==False)&(CentralDifference[1]==True)):
            dEx_dx=((((IncBeam[4,linePositionAbove_IncBeam]+IntField[4,linePositionAbove_IntField])-(IncBeam[4,linePositionCurrent_IncBeam]+IntField[4,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[5,linePositionAbove_IncBeam]+IntField[5,linePositionAbove_IntField])-(IncBeam[5,linePositionCurrent_IncBeam]+IntField[5,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEy_dx=((((IncBeam[6,linePositionAbove_IncBeam]+IntField[6,linePositionAbove_IntField])-(IncBeam[6,linePositionCurrent_IncBeam]+IntField[6,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[7,linePositionAbove_IncBeam]+IntField[7,linePositionAbove_IntField])-(IncBeam[7,linePositionCurrent_IncBeam]+IntField[7,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEz_dx=((((IncBeam[8,linePositionAbove_IncBeam]+IntField[8,linePositionAbove_IntField])-(IncBeam[8,linePositionCurrent_IncBeam]+IntField[8,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[9,linePositionAbove_IncBeam]+IntField[9,linePositionAbove_IntField])-(IncBeam[9,linePositionCurrent_IncBeam]+IntField[9,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            Force[i,4]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dx)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dx)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dx))
        if ((CentralDifference[0]==True)&(CentralDifference[1]==False)):
            dEx_dx=((((IncBeam[4,linePositionCurrent_IncBeam]+IntField[4,linePositionCurrent_IntField])-(IncBeam[4,linePositionBelow_IncBeam]+IntField[4,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[5,linePositionCurrent_IncBeam]+IntField[5,linePositionCurrent_IntField])-(IncBeam[5,linePositionBelow_IncBeam]+IntField[5,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEy_dx=((((IncBeam[6,linePositionCurrent_IncBeam]+IntField[6,linePositionCurrent_IntField])-(IncBeam[6,linePositionBelow_IncBeam]+IntField[6,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[7,linePositionCurrent_IncBeam]+IntField[7,linePositionCurrent_IntField])-(IncBeam[7,linePositionBelow_IncBeam]+IntField[7,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEz_dx=((((IncBeam[8,linePositionCurrent_IncBeam]+IntField[8,linePositionCurrent_IntField])-(IncBeam[8,linePositionBelow_IncBeam]+IntField[8,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[9,linePositionCurrent_IncBeam]+IntField[9,linePositionCurrent_IntField])-(IncBeam[9,linePositionBelow_IncBeam]+IntField[9,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            Force[i,4]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dx)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dx)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dx))
            
        #Calculate the force in y direction
        linePositionBelow_IncBeam=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((3*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]-((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionCurrent_IncBeam=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionAbove_IncBeam=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]+((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((3*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionBelow_IntField=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((3*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]-((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionCurrent_IntField=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionAbove_IntField=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]+((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((3*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        CentralDifference=np.ones([2],dtype=bool) #To check if lower and upper part of Central difference is available
        if((len(linePositionBelow_IncBeam)==0)or(len(linePositionBelow_IntField)==0)): #Is there a value below?
            CentralDifference[0]=False
        if((len(linePositionAbove_IncBeam)==0)or(len(linePositionAbove_IntField)==0)): #Is there a value above?
            CentralDifference[1]=False
        if ((CentralDifference[0]==True)&(CentralDifference[1]==True)):
            dEx_dy=((((IncBeam[4,linePositionAbove_IncBeam]+IntField[4,linePositionAbove_IntField])-(IncBeam[4,linePositionBelow_IncBeam]+IntField[4,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[5,linePositionAbove_IncBeam]+IntField[5,linePositionAbove_IntField])-(IncBeam[5,linePositionBelow_IncBeam]+IntField[5,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            dEy_dy=((((IncBeam[6,linePositionAbove_IncBeam]+IntField[6,linePositionAbove_IntField])-(IncBeam[6,linePositionBelow_IncBeam]+IntField[6,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[7,linePositionAbove_IncBeam]+IntField[7,linePositionAbove_IntField])-(IncBeam[7,linePositionBelow_IncBeam]+IntField[7,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            dEz_dy=((((IncBeam[8,linePositionAbove_IncBeam]+IntField[8,linePositionAbove_IntField])-(IncBeam[8,linePositionBelow_IncBeam]+IntField[8,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[9,linePositionAbove_IncBeam]+IntField[9,linePositionAbove_IntField])-(IncBeam[9,linePositionBelow_IncBeam]+IntField[9,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            Force[i,5]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dy)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dy)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dy))
        if ((CentralDifference[0]==False)&(CentralDifference[1]==True)):
            dEx_dy=((((IncBeam[4,linePositionAbove_IncBeam]+IntField[4,linePositionAbove_IntField])-(IncBeam[4,linePositionCurrent_IncBeam]+IntField[4,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[5,linePositionAbove_IncBeam]+IntField[5,linePositionAbove_IntField])-(IncBeam[5,linePositionCurrent_IncBeam]+IntField[5,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEy_dy=((((IncBeam[6,linePositionAbove_IncBeam]+IntField[6,linePositionAbove_IntField])-(IncBeam[6,linePositionCurrent_IncBeam]+IntField[6,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[7,linePositionAbove_IncBeam]+IntField[7,linePositionAbove_IntField])-(IncBeam[7,linePositionCurrent_IncBeam]+IntField[7,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEz_dy=((((IncBeam[8,linePositionAbove_IncBeam]+IntField[8,linePositionAbove_IntField])-(IncBeam[8,linePositionCurrent_IncBeam]+IntField[8,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[9,linePositionAbove_IncBeam]+IntField[9,linePositionAbove_IntField])-(IncBeam[9,linePositionCurrent_IncBeam]+IntField[9,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            Force[i,5]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dy)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dy)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dy))
        if ((CentralDifference[0]==True)&(CentralDifference[1]==False)):
            dEx_dy=((((IncBeam[4,linePositionCurrent_IncBeam]+IntField[4,linePositionCurrent_IntField])-(IncBeam[4,linePositionBelow_IncBeam]+IntField[4,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[5,linePositionCurrent_IncBeam]+IntField[5,linePositionCurrent_IntField])-(IncBeam[5,linePositionBelow_IncBeam]+IntField[5,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEy_dy=((((IncBeam[6,linePositionCurrent_IncBeam]+IntField[6,linePositionCurrent_IntField])-(IncBeam[6,linePositionBelow_IncBeam]+IntField[6,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[7,linePositionCurrent_IncBeam]+IntField[7,linePositionCurrent_IntField])-(IncBeam[7,linePositionBelow_IncBeam]+IntField[7,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEz_dy=((((IncBeam[8,linePositionCurrent_IncBeam]+IntField[8,linePositionCurrent_IntField])-(IncBeam[8,linePositionBelow_IncBeam]+IntField[8,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[9,linePositionCurrent_IncBeam]+IntField[9,linePositionCurrent_IntField])-(IncBeam[9,linePositionBelow_IncBeam]+IntField[9,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            Force[i,5]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dy)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dy)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dy))
        
        """#Calculate the force in z direction
        linePositionBelow_IncBeam=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((3*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]-((1*DipoleSep)/2))))
        linePositionCurrent_IncBeam=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionAbove_IncBeam=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]+((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((3*DipoleSep)/2))))
        linePositionBelow_IntField=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((3*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]-((1*DipoleSep)/2))))
        linePositionCurrent_IntField=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
        linePositionAbove_IntField=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]+((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((3*DipoleSep)/2))))
        CentralDifference=np.ones([2],dtype=bool) #To check if lower and upper part of Central difference is available
        if((len(linePositionBelow_IncBeam)==0)or(len(linePositionBelow_IntField)==0)): #Is there a value below?
            CentralDifference[0]=False
        if((len(linePositionAbove_IncBeam)==0)or(len(linePositionAbove_IntField)==0)): #Is there a value above?
            CentralDifference[1]=False
        if ((CentralDifference[0]==True)&(CentralDifference[1]==True)):
            dEx_dz=((((IncBeam[4,linePositionAbove_IncBeam]+IntField[4,linePositionAbove_IntField])-(IncBeam[4,linePositionBelow_IncBeam]+IntField[4,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[5,linePositionAbove_IncBeam]+IntField[5,linePositionAbove_IntField])-(IncBeam[5,linePositionBelow_IncBeam]+IntField[5,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            dEy_dz=((((IncBeam[6,linePositionAbove_IncBeam]+IntField[6,linePositionAbove_IntField])-(IncBeam[6,linePositionBelow_IncBeam]+IntField[6,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[7,linePositionAbove_IncBeam]+IntField[7,linePositionAbove_IntField])-(IncBeam[7,linePositionBelow_IncBeam]+IntField[7,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            dEz_dz=((((IncBeam[8,linePositionAbove_IncBeam]+IntField[8,linePositionAbove_IntField])-(IncBeam[8,linePositionBelow_IncBeam]+IntField[8,linePositionBelow_IntField]))/(2*DipoleSep))-(1j*(((IncBeam[9,linePositionAbove_IncBeam]+IntField[9,linePositionAbove_IntField])-(IncBeam[9,linePositionBelow_IncBeam]+IntField[9,linePositionBelow_IntField]))/(2*DipoleSep)))) #Complex conjugate already implemented
            Force[i,6]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dz)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dz)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dz))
        if ((CentralDifference[0]==False)&(CentralDifference[1]==True)):
            dEx_dz=((((IncBeam[4,linePositionAbove_IncBeam]+IntField[4,linePositionAbove_IntField])-(IncBeam[4,linePositionCurrent_IncBeam]+IntField[4,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[5,linePositionAbove_IncBeam]+IntField[5,linePositionAbove_IntField])-(IncBeam[5,linePositionCurrent_IncBeam]+IntField[5,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEy_dz=((((IncBeam[6,linePositionAbove_IncBeam]+IntField[6,linePositionAbove_IntField])-(IncBeam[6,linePositionCurrent_IncBeam]+IntField[6,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[7,linePositionAbove_IncBeam]+IntField[7,linePositionAbove_IntField])-(IncBeam[7,linePositionCurrent_IncBeam]+IntField[7,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEz_dz=((((IncBeam[8,linePositionAbove_IncBeam]+IntField[8,linePositionAbove_IntField])-(IncBeam[8,linePositionCurrent_IncBeam]+IntField[8,linePositionCurrent_IntField]))/(DipoleSep))-(1j*(((IncBeam[9,linePositionAbove_IncBeam]+IntField[9,linePositionAbove_IntField])-(IncBeam[9,linePositionCurrent_IncBeam]+IntField[9,linePositionCurrent_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            Force[i,6]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dz)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dz)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dz))
        if ((CentralDifference[0]==True)&(CentralDifference[1]==False)):
            dEx_dz=((((IncBeam[4,linePositionCurrent_IncBeam]+IntField[4,linePositionCurrent_IntField])-(IncBeam[4,linePositionBelow_IncBeam]+IntField[4,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[5,linePositionCurrent_IncBeam]+IntField[5,linePositionCurrent_IntField])-(IncBeam[5,linePositionBelow_IncBeam]+IntField[5,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEy_dz=((((IncBeam[6,linePositionCurrent_IncBeam]+IntField[6,linePositionCurrent_IntField])-(IncBeam[6,linePositionBelow_IncBeam]+IntField[6,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[7,linePositionCurrent_IncBeam]+IntField[7,linePositionCurrent_IntField])-(IncBeam[7,linePositionBelow_IncBeam]+IntField[7,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            dEz_dz=((((IncBeam[8,linePositionCurrent_IncBeam]+IntField[8,linePositionCurrent_IntField])-(IncBeam[8,linePositionBelow_IncBeam]+IntField[8,linePositionBelow_IntField]))/(DipoleSep))-(1j*(((IncBeam[9,linePositionCurrent_IncBeam]+IntField[9,linePositionCurrent_IntField])-(IncBeam[9,linePositionBelow_IncBeam]+IntField[9,linePositionBelow_IntField]))/(DipoleSep)))) #Complex conjugate already implemented
            Force[i,6]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dz)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dz)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dz))"""
        Force[i,6] = 0
            
        Force[i,3]=pow((pow(Force[i,4],2)+pow(Force[i,5],2)+pow(Force[i,6],2)),0.5)
    
    return Force

#Preliminary Variables	
CorrectionFactor=0.0838
BeamWidth=0.6 #In micro m
Temperature=20 #Degrees C
Power=5e-3 #In Watts
MediumDielectricConstant=87.740-(0.40008*Temperature)+(9.398e-4*(Temperature**2))-(1.410e-6*(Temperature**3))
ElectricFieldStrength=ElectricFieldStrengthCalc(MediumDielectricConstant,Power,BeamWidth) #V/m
nu = 8.891e-4
r = 1e-6

#Preliminary Dynamic variables
t_0 = 0
t_end = 1
t_step = 5e-4

#Start point of beam 1
x_beam1, y_beam1, z_beam1 = 0,0,0

#Start point of beam 2
x_beam2, y_beam2, z_beam2 = 0,2,0

#Particle start point start point
x_position, y_position, z_position = 0,0,0

x_beam1 += -x_position
y_beam1 += -y_position
z_beam1 += -z_position
x_beam2 += -x_position
y_beam2 += -y_position
z_beam2 += -z_position

#Array to track particle position
PPositionArray = np.zeros([1,5])
StartTime=time.clock()
while t_0 <= t_end:

    #Perform the DDA Calculations and calculate forces

    DipPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'DipPol-Y')  
    IntFPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IntField-Y')
    BeamPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IncBeam-Y')
    ForcePathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'RadForce-Y')
    print('Processing time: '+str(t_0))
    
    #Calculate incident beam due to first beam
    callString=".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl 15 -m 1.18339034696 0 -lambda 1.064 -prop 0 0 1 -beam barton5 "+str(BeamWidth)+" "+str(x_beam1)+" "+str(y_beam1)+" "+str(z_beam1)+" -store_beam" #The script for performing the DDA calculations
    subprocess.call(callString,shell=True)
    Beam1File = sorted(glob.glob(BeamPathInput))[-1]
    Beam1Raw=np.transpose(np.loadtxt(Beam1File, skiprows=1))
    try:
        shutil.rmtree(Beam1File.replace(os.sep+'IncBeam-Y',''))
    except:
        print('Cannot Delete')
        
    #Calculate incident beam due to second beam
    callString=".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl 15 -m 1.18339034696 0 -lambda 1.064 -prop 0 0 1 -beam barton5 "+str(BeamWidth)+" "+str(x_beam2)+" "+str(y_beam2)+" "+str(z_beam2)+" -store_beam" #The script for performing the DDA calculations
    subprocess.call(callString,shell=True)
    Beam2File = sorted(glob.glob(BeamPathInput))[-1]
    Beam2Raw=np.transpose(np.loadtxt(Beam2File, skiprows=1))
    try:
        shutil.rmtree(Beam2File.replace(os.sep+'IncBeam-Y',''))
    except:
        print('Cannot Delete')
        
    #Calculate the new beam file
    IncBeamRaw=np.zeros([np.size(Beam1Raw,axis=0),np.size(Beam1Raw,axis=1)])
    for i in range(0,3): #Bring across the x,y,z positions
        IncBeamRaw[i,:]=Beam1Raw[i,:]
    for i in range(4,10): #Bring across the total E_inc beam
        IncBeamRaw[i,:]=np.add(Beam1Raw[i,:],Beam2Raw[i,:])
        
    IncBeamRaw[3,:]=np.add((np.add((np.power((np.abs(np.add(IncBeamRaw[4,:],(1j*IncBeamRaw[5,:])))),2)),(np.power((np.abs(np.add(IncBeamRaw[6,:],(1j*IncBeamRaw[7,:])))),2)))),(np.power((np.abs(np.add(IncBeamRaw[8,:],(1j*IncBeamRaw[9,:])))),2)))
    DualBeamPath = str(os.getcwd())+str(os.sep+'DualBeam')	
    with open(DualBeamPath, 'wb') as f:
        f.write(b'x y z |Einc|^2 Eincx.r Eincx.i Eincy.r Eincy.i Eincz.r Eincz.i\n')
        np.savetxt(f, np.transpose(IncBeamRaw), fmt='%.11f', delimiter=' ')
    
    
    callString=".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl 15 -m 1.18339034696 0 -lambda 1.064 -prop 0 0 1 -sym enf -beam read DualBeam -store_dip_pol -store_int_field" #The script for performing the DDA calculations
    print(".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl 15 -m 1.18339034696 0 -lambda 1.064 -prop 0 0 1 -sym enf -beam read DualBeam -store_dip_pol -store_int_field")
    subprocess.call(callString,shell=True)
    DipFiles, IntFFiles= sorted(glob.glob(DipPathInput))[-1], sorted(glob.glob(IntFPathInput))[-1] #File containing the paths to each DipPol, IntField file
    DipPolRaw=np.transpose(np.loadtxt(DipFiles, skiprows=1))
    IntFieldRaw=np.transpose(np.loadtxt(IntFFiles, skiprows=1))
    DipoleSeperation=DipSep(DipPolRaw[0,:])
    StartTime_OurCalc=time.clock()
    CalculatedForce=Forces(DipPolRaw,IntFieldRaw,IncBeamRaw,DipoleSeperation)
    EndTime_OurCalc=time.clock()
    
    
    '''#SAVE CALCULATED FORCES
    FFiles="Force"
    with open(FFiles,'wb') as f:
        f.write(b'x y z |F|^2 Fx Fy Fz \n')
        np.savetxt(f,CalculatedForce, fmt='%e',delimiter=' ')'''
    
    #This section is where we look at the ADDA Calculated Forces
    EstimatedParticleForce1=np.array([[np.sum(CalculatedForce[:,4])],[np.sum(CalculatedForce[:,5])],[np.sum(CalculatedForce[:,6])]])
    EstimatedParticleForce2=OurForceConversion(EstimatedParticleForce1,CorrectionFactor,ElectricFieldStrength) #Convert to ADDA and SI

    #Generate the Brownian "Force"
    Drag_Coefficient = DragCoef(nu,r)
    B_x, B_y, B_z = BrownianForce(Drag_Coefficient, Temperature), BrownianForce(Drag_Coefficient, Temperature), 0 #BrownianForce(Drag_Coefficient, Temperature)

    EstimatedParticleForce3 = np.array([[EstimatedParticleForce2[0]+B_x],[EstimatedParticleForce2[1]+B_y],[EstimatedParticleForce2[2]+B_z]])
                    
    #Calculate the change in position of the particle
    x = PositionChange(EstimatedParticleForce3[0], Drag_Coefficient, t_step)
    y = PositionChange(EstimatedParticleForce3[1], Drag_Coefficient, t_step)
    z = PositionChange(EstimatedParticleForce3[2], Drag_Coefficient, t_step)
    x_beam1 += -x[0,0]
    y_beam1 += -y[0,0]
    #z_beam1 += -z[0,0]
    x_beam2 += -x[0,0]
    y_beam2 += -y[0,0]
    #z_beam2 += -z[0,0]
    x_position +=x[0,0]
    y_position +=y[0,0]
    #z_position +=z[0,0]
    t_0 += t_step
    PPositionArray = np.append(PPositionArray, np.array([[t_0,x_position,y_position,EstimatedParticleForce3[0],EstimatedParticleForce3[1]]]),axis=0)
                       
    #Use to delete the files after processing
    try:
        shutil.rmtree(DipFiles.replace(os.sep+'DipPol-Y',''))
    except:
        print('Cannot Delete')

    EndTime=time.clock()
    TimeRecordings=np.array([[(EndTime-StartTime)]])        
    np.savetxt('ParticlePositions', PPositionArray, fmt='%e', delimiter=' ')
    TimeLogPath = str(os.getcwd())+str(os.sep+'TimeLog')	
    with open(TimeLogPath, 'wb') as f:
        f.write(b'Time(s)\n')
        np.savetxt(f, TimeRecordings, fmt='%.10f', delimiter=' ')

