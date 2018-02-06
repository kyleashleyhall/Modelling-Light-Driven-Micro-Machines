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
Initial_dpl=19
Final_dpl=20
Step_dpl=1

Initial_lambdaValue=0.9
Final_lambdaValue=1.1
Step_lambdaValue=0.1

Initial_refractiveIndexValue=1.4
Final_refractiveIndexValue=1.65
Step_refractiveIndexValue=0.05

Initial_particleDiameterValue=1
Final_particleDiameterValue=2.1
Step_particleDiameterValue=0.1


#Perform the DDA Calculations and calculate forces
DipPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'DipPol-Y')  
IntFPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IntField-Y')
BeamPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IncBeam-Y')
ForcePathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'RadForce-Y')
dpl=Initial_dpl
CorrectionFactors=np.zeros([462,5]) #Edit me if varying Preliminary Variables
CalculationTimes=np.zeros([(Final_dpl-Initial_dpl),2])
iterator=0
while (dpl<Final_dpl):
    StartTime=time.clock()
    lambdaValue=Initial_lambdaValue
    while(lambdaValue<Final_lambdaValue):
        refractiveIndexValue=Initial_refractiveIndexValue
        while(refractiveIndexValue<Final_refractiveIndexValue):
            particleDiameterValue=Initial_particleDiameterValue
            while(particleDiameterValue<Final_particleDiameterValue):
                
                callString=".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size "+str(particleDiameterValue)+" -dpl "+str(dpl)+" -lambda "+str(lambdaValue)+" -prop 0 0 1 -m "+str(refractiveIndexValue)+" 0 -store_beam -store_dip_pol -store_int_field -store_force" #The script for performing the DDA calculations
                subprocess.call(callString,shell=True)
                DipFiles, IntFFiles, BeamFiles, ForceFiles = sorted(glob.glob(DipPathInput))[-1], sorted(glob.glob(IntFPathInput))[-1], sorted(glob.glob(BeamPathInput))[-1], sorted(glob.glob(ForcePathInput))[-1] #File containing the paths to each DipPol, IntField file
                FFiles = DipFiles.replace('DipPol-Y','CalculatedForces')
                DipPolRaw=np.transpose(np.loadtxt(DipFiles, skiprows=1))
                IntFieldRaw=np.transpose(np.loadtxt(IntFFiles, skiprows=1))
                IncBeamRaw=np.transpose(np.loadtxt(BeamFiles, skiprows=1))
                DipoleSeperation=DipSep(DipPolRaw[0,:])
                CalculatedForce=Forces(DipPolRaw,IntFieldRaw,IncBeamRaw,DipoleSeperation)
    
    
                #SAVE CALCULATED FORCES
                with open(FFiles,'wb') as f:
                    f.write(b'x y z |F|^2 Fx Fy Fz \n')
                    np.savetxt(f,CalculatedForce, fmt='%e',delimiter=' ')
    
                #This section is where we look at the ADDA Calculated Forces
                EstimatedParticleForce=np.array([[np.sum(CalculatedForce[:,4])],[np.sum(CalculatedForce[:,5])],[np.sum(CalculatedForce[:,6])]])
                ADDADipoleForceFile = np.loadtxt(ForceFiles, skiprows=1) #Load the ADDA Dipole Forces File
                ADDAParticleForce = np.array([[np.sum(ADDADipoleForceFile[:,4])],[np.sum(ADDADipoleForceFile[:,5])],[np.sum(ADDADipoleForceFile[:,6])]]) #Save the ADDA Particle Forces to memory
                CorrectionFactors[iterator,0] = dpl
                CorrectionFactors[iterator,1] = lambdaValue
                CorrectionFactors[iterator,2] = refractiveIndexValue
                CorrectionFactors[iterator,3] = particleDiameterValue
                CorrectionFactors[iterator,4] = ADDAParticleForce[2]/EstimatedParticleForce[2]
                iterator+=1
                particleDiameterValue+=Step_particleDiameterValue
                
                #Use to delete the files after processing
                try:
                    shutil.rmtree(FFiles.replace(os.sep+'CalculatedForces',''))
                except:
                    print('Cannot Delete')
                
            refractiveIndexValue+=Step_refractiveIndexValue
            
        lambdaValue+=Step_lambdaValue
        
    #Include in dpl script
    EndTime=time.clock()
    CalculationTimes[(dpl-Initial_dpl),0] = dpl
    CalculationTimes[(dpl-Initial_dpl),1] = EndTime-StartTime
    dpl=dpl+Step_dpl


CorrectionFactorsPath = str(os.getcwd())+str(os.sep+'CorrectionFactors') #|(F_ADDA - F_Calc)|/F_ADDA	
with open(CorrectionFactorsPath, 'wb') as f:
    f.write(b'dpl Lambda(micro m) RefractiveIndex ParticleDiameter(micro m) CorrectionFactor\n')
    np.savetxt(f, CorrectionFactors, fmt='%.10f', delimiter=' ')
TimeLogPath = str(os.getcwd())+str(os.sep+'TimeLog')	
with open(TimeLogPath, 'wb') as f:
    f.write(b'dpl Time\n')
    np.savetxt(f, CalculationTimes, fmt='%.10f', delimiter=' ')