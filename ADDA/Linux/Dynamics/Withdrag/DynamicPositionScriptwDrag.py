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
import scipy.constants as constants

def deltapos(F,m,t,v):
    '''Solution to F=ma, deltapos stands for change in position'''
    return float(F*(t**2)/(2*m) + v*t)
    
def deltavel(F,m,t):
    return F*t/m

def DragForce(nu, r, vel):
    return 6*np.pi*nu*r*vel
        
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
        
        #Calculate the force in z direction
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
            Force[i,6]=0.5*np.real(((DipPol[4,i]+(1j*DipPol[5,i]))*dEx_dz)+((DipPol[6,i]+(1j*DipPol[7,i]))*dEy_dz)+((DipPol[8,i]+(1j*DipPol[9,i]))*dEz_dz))
            
        Force[i,3]=pow((pow(Force[i,4],2)+pow(Force[i,5],2)+pow(Force[i,6],2)),0.5)
    
    return Force

#Preliminary Variables	
ElectricFieldStrength=1 #V/m
CorrectionFactor=0.0734357589	

m = 4.2e-12 #Polystyrene bead density*1 micrometer radius
nu = 8.891e-4
r = 1e-6 

#Preliminary Dynamic variables
t_0 = 0
t_end = 0.4
t_step = 0.1
x_beam, y_beam, z_beam = 0,0,0
v_x, v_y, v_z = 0,0,0

#Array to track particle position
N = (t_end - t_0) / t_step
PPositionArray = np.zeros([1,4])
while t_0 <= t_end:

    #Perform the DDA Calculations and calculate forces

    DipPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'DipPol-Y')  
    IntFPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IntField-Y')
    BeamPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IncBeam-Y')
    ForcePathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'RadForce-Y')
    #TimeRecordings=np.zeros([(Final_dpl-Initial_dpl),3])
    #CalculationTimes=np.zeros([(Final_dpl-Initial_dpl),2])
    print('Processing time: '+str(t_0))
    callString=".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl 15 -lambda 1 -prop 0 0 1 -beam barton5 1 "+str(x_beam)+" "+str(y_beam)+" "+str(z_beam)+" -store_beam -store_dip_pol -store_int_field" #The script for performing the DDA calculations
    print(".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl 15 -lambda 1 -prop 0 0 1 -beam barton5 1 "+str(x_beam)+" "+str(y_beam)+" "+str(z_beam)+" -store_beam -store_dip_pol -store_int_field")
    StartTime_ADDA=time.clock()
    subprocess.call(callString,shell=True)
    EndTime_ADDA=time.clock()
    DipFiles, IntFFiles, BeamFiles = sorted(glob.glob(DipPathInput))[-1], sorted(glob.glob(IntFPathInput))[-1], sorted(glob.glob(BeamPathInput))[-1] #File containing the paths to each DipPol, IntField file
    DipPolRaw=np.transpose(np.loadtxt(DipFiles, skiprows=1))
    IntFieldRaw=np.transpose(np.loadtxt(IntFFiles, skiprows=1))
    IncBeamRaw=np.transpose(np.loadtxt(BeamFiles, skiprows=1))
    DipoleSeperation=DipSep(DipPolRaw[0,:])
    StartTime_OurCalc=time.clock()
    CalculatedForce=Forces(DipPolRaw,IntFieldRaw,IncBeamRaw,DipoleSeperation)
    EndTime_OurCalc=time.clock()
    
    
    #SAVE CALCULATED FORCES
#    with open(FFiles,'wb') as f:
#        f.write(b'x y z |F|^2 Fx Fy Fz \n')
#        np.savetxt(f,CalculatedForce, fmt='%e',delimiter=' ')
    
    #This section is where we look at the ADDA Calculated Forces
    EstimatedParticleForce1=np.array([[np.sum(CalculatedForce[:,4])],[np.sum(CalculatedForce[:,5])],[np.sum(CalculatedForce[:,6])]])
    EstimatedParticleForce2=OurForceConversion(EstimatedParticleForce1,CorrectionFactor,ElectricFieldStrength) #Convert to ADDA and SI
    v_x += deltavel(EstimatedParticleForce2[0],m,t_step)
    v_y += deltavel(EstimatedParticleForce2[1],m,t_step)
    v_z += deltavel(EstimatedParticleForce2[2],m,t_step)
    print(v_x,v_y,v_z)
    F_dragx = DragForce(nu, r, v_x)
    F_dragy = DragForce(nu, r, v_y)
    F_dragz = DragForce(nu, r, v_z)
    print(F_dragx,F_dragy,F_dragz)
    EstimatedParticleForce3 = np.array([[EstimatedParticleForce2[0]-F_dragx],[EstimatedParticleForce2[1]-F_dragy],[EstimatedParticleForce2[2]-F_dragz]])
#==============================================================================
#         ADDADipoleForceFile = np.loadtxt(ForceFiles, skiprows=1) #Load the ADDA Dipole Forces File
#         ADDAParticleForce = np.array([[np.sum(ADDADipoleForceFile[:,4])],[np.sum(ADDADipoleForceFile[:,5])],[np.sum(ADDADipoleForceFile[:,6])]]) #Save the ADDA Particle Forces to memory
#         ForceError[(dpl-Initial_dpl),0] = dpl 
#         ForceError[(dpl-Initial_dpl),1] = pow((pow(ADDAParticleForce[0]-EstimatedParticleForce[0],2)),0.5)
#         ForceError[(dpl-Initial_dpl),2] = pow((pow(ADDAParticleForce[1]-EstimatedParticleForce[1],2)),0.5)
#         ForceError[(dpl-Initial_dpl),3] = pow((pow(ADDAParticleForce[2]-EstimatedParticleForce[2],2)),0.5)
#         ForceError[(dpl-Initial_dpl),4] = pow((pow(ADDAParticleForce[0]-EstimatedParticleForce[0],2)),0.5)/ADDAParticleForce[0]
#         ForceError[(dpl-Initial_dpl),5] = pow((pow(ADDAParticleForce[1]-EstimatedParticleForce[1],2)),0.5)/ADDAParticleForce[1]
#         ForceError[(dpl-Initial_dpl),6] = pow((pow(ADDAParticleForce[2]-EstimatedParticleForce[2],2)),0.5)/ADDAParticleForce[2]
#         TimeRecordings[(dpl-Initial_dpl),0] = dpl
#         TimeRecordings[(dpl-Initial_dpl),1] = EndTime_ADDA-StartTime_ADDA
#         TimeRecordings[(dpl-Initial_dpl),2] = EndTime_OurCalc-StartTime_OurCalc
#==============================================================================
                    
    #Calculate the change in position of the particle
    x = deltapos(EstimatedParticleForce3[0],m,t_step, v_x)
    y = deltapos(EstimatedParticleForce3[1],m,t_step, v_y)
    z = deltapos(EstimatedParticleForce3[2],m,t_step, v_z)
    x_beam += -x
    y_beam += -y
    z_beam += -z
    t_0 += t_step
    PPositionArray = np.append(PPositionArray, np.array([[t_0,-x_beam,-y_beam,-z_beam]]),axis=0)
                       
    #Use to delete the files after processing
    try:
        shutil.rmtree(DipFiles.replace(os.sep+'DipPol-Y',''))
    except:
        print('Cannot Delete')
        
np.savetxt('ParticlePositions', PPositionArray, fmt='%e', delimiter=' ')
#==============================================================================
# ForceErrorPath = str(os.getcwd())+str(os.sep+'ForceError') #|(F_ADDA - F_Calc)|/F_ADDA	
# with open(ForceErrorPath, 'wb') as f:
#     f.write(b'dpl |F_ADDA(x)-F_Calc(x)| |F_ADDA(y)-F_Calc(y)| |F_ADDA(z)-F_Calc(z)| |F_ADDA(x)-F_Calc(x)|/F_ADDA(x) |F_ADDA(y)-F_Calc(y)|/F_ADDA(y) |F_ADDA(z)-F_Calc(z)|/F_ADDA(z)\n')
#     np.savetxt(f, ForceError, fmt='%.10f', delimiter=' ')
# TimeLogPath = str(os.getcwd())+str(os.sep+'TimeLog')	
# with open(TimeLogPath, 'wb') as f:
#     f.write(b'dpl ADDATime OurCalcTime\n')
#     np.savetxt(f, TimeRecordings, fmt='%.10f', delimiter=' ')
#==============================================================================
