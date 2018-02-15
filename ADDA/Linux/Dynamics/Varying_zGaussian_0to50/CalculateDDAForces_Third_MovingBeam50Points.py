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

def ADDAForceConversion(Force,ElectricFieldStrength):

	Force=Force*(((ElectricFieldStrength)**2)/((constants.c)**2))*(10**(-5))
	
	return(Force)
	
def OurForceConversion(Force,CorrectionFactor,ElectricFieldStrength):
    
   return(ADDAForceConversion((Force*CorrectionFactor),ElectricFieldStrength))

def DipSep(Singleaxis):
    dx = np.zeros([len(Singleaxis)-1])
    for i in range(len(Singleaxis)-1):
        dx[i] = Singleaxis[i+1] - Singleaxis[i]
        if (dx[i]==0):
            dx[i]=np.inf
    return min(np.absolute(dx))
    
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
Initial_dpl=15
Final_dpl=16
Step_dpl=1
x=0
y=0
z=0
Initial_z=0
Step_z=0.5
Final_z=50

ElectricFieldStrength=1e8 #V/m --We would rather the force be too large than too small for this test - Marius van Laar 12/02/2018 16:45
CorrectionFactor=0.0734357589

z=Initial_z #Set the value of y to the initial value

#Perform the DDA Calculations and calculate forces
DipPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'DipPol-Y')  
IntFPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IntField-Y')
BeamPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IncBeam-Y')
dpl=Initial_dpl
TimeRecordings=np.zeros([(((Final_z-Initial_z)//Step_dpl)+1),5])
while (dpl<Final_dpl):
    while (z<Final_z):
        print('Processing dpl: '+str(dpl))
        callString=".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl "+str(dpl)+" -lambda 0.9 -m 1.1859519224 0 -prop 0 0 1 -beam barton5 1 "+str(x)+" "+str(y)+" "+str(-z)+" -store_beam -store_dip_pol -store_int_field" #The script for performing the DDA calculations
        print(".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl "+str(dpl)+" -lambda 0.9 -m 1.1859519224 0 -prop 0 0 1 -beam barton5 1 "+str(x)+" "+str(y)+" "+str(-z)+" -store_beam -store_dip_pol -store_int_field")
        StartTime_ADDA=time.clock()
        subprocess.call(callString,shell=True)
        EndTime_ADDA=time.clock()
        DipFiles, IntFFiles, BeamFiles = sorted(glob.glob(DipPathInput))[-1], sorted(glob.glob(IntFPathInput))[-1], sorted(glob.glob(BeamPathInput))[-1] #File containing the paths to each DipPol, IntField file
        FFiles = DipFiles.replace('DipPol-Y','CalculatedForces')
        DipPolRaw=np.transpose(np.loadtxt(DipFiles, skiprows=1))
        IntFieldRaw=np.transpose(np.loadtxt(IntFFiles, skiprows=1))
        IncBeamRaw=np.transpose(np.loadtxt(BeamFiles, skiprows=1))
        DipoleSeperation=DipSep(DipPolRaw[0,:])
        StartTime_OurCalc=time.clock()
        CalculatedForce=Forces(DipPolRaw,IntFieldRaw,IncBeamRaw,DipoleSeperation)
        Force_x=OurForceConversion(np.sum(CalculatedForce[:,4]),CorrectionFactor,ElectricFieldStrength)
        Force_y=OurForceConversion(np.sum(CalculatedForce[:,5]),CorrectionFactor,ElectricFieldStrength)
        Force_z=OurForceConversion(np.sum(CalculatedForce[:,6]),CorrectionFactor,ElectricFieldStrength)
        EndTime_OurCalc=time.clock()
    
    
        #SAVE CALCULATED FORCES
        with open(FFiles,'wb') as f:
            f.write(b'x y z |F|^2 Fx Fy Fz \n')
            np.savetxt(f,CalculatedForce, fmt='%e',delimiter=' ')
    
        #This section is where we look at the ADDA Calculated Forces
        EstimatedParticleForce=np.array([[x],[y],[z],[Force_x],[Force_y],[Force_z]])
        try:
            AllForces = np.hstack([AllForces,EstimatedParticleForce])
        except:
            AllForces = EstimatedParticleForce
    
        TimeRecordings[(((z-Initial_z)//Step_dpl)),0] = x
        TimeRecordings[(((z-Initial_z)//Step_dpl)),1] = y
        TimeRecordings[(((z-Initial_z)//Step_dpl)),2] = z
        TimeRecordings[(((z-Initial_z)//Step_dpl)),3] = EndTime_ADDA-StartTime_ADDA
        TimeRecordings[(((z-Initial_z)//Step_dpl)),4] = EndTime_OurCalc-StartTime_OurCalc
    
        #Use to delete the files after processing
        '''try:
            shutil.rmtree(FFiles.replace(os.sep+'CalculatedForces',''))
        except:
            print('Cannot Delete')'''
        z=z+Step_z
    dpl=dpl+Step_dpl

TimeLogPath = str(os.getcwd())+str(os.sep+'TimeLog')	
with open(TimeLogPath, 'wb') as f:
    f.write(b'x y z ADDATime OurCalcTime\n')
    np.savetxt(f, TimeRecordings, fmt='%.10f', delimiter=' ')
AllForcesPath = str(os.getcwd())+str(os.sep+'TotalForce')
np.savetxt(AllForcesPath,np.transpose(AllForces), fmt='%.10f',delimiter=' ')    