# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:06:32 2018

@author: Kyle
"""

import numpy as np


def RotateBeam(BeamFocusPoint,BeamProp,CentreOfMass,phi,theta,psi):
    
    #Construct the rotation matrix
    a_11=((np.cos(psi))*(np.cos(phi)))-((np.cos(theta))*(np.sin(phi))*(np.sin(psi)))
    a_12=((np.cos(psi))*(np.sin(phi)))+((np.cos(theta))*(np.cos(phi))*(np.sin(psi)))
    a_13=(np.sin(psi))*(np.sin(theta))
    a_21=-((np.sin(psi))*(np.cos(phi)))-((np.cos(theta))*(np.sin(phi))*(np.cos(psi)))
    a_22=-((np.sin(psi))*(np.sin(phi)))+((np.cos(theta))*(np.cos(phi))*(np.cos(psi)))
    a_23=(np.cos(psi))*(np.sin(theta))
    a_31=(np.sin(theta))*(np.sin(phi))
    a_32=-((np.sin(theta))*(np.cos(phi)))
    a_33=np.cos(theta)
    RotationMatrix_A=np.array([[a_11,a_12,a_13],[a_21,a_22,a_23],[a_31,a_32,a_33]])
    
    CoMTransposed_BeamFocusPoint=BeamFocusPoint-CentreOfMass #Move the points so the Centre of Mass is at (0,0,0)
    
    CoMTransposedPropPoint=CoMTransposed_BeamFocusPoint+BeamProp #Calculate the point above so that we can rotate the propogation direction
    
    CoMTransposed_Rotated_BeamFocusPoint=np.dot(RotationMatrix_A,CoMTransposed_BeamFocusPoint) #Calculate the rotated focus point
    
    CoMTransposed_Rotated_PropPoint=np.dot(RotationMatrix_A,CoMTransposedPropPoint) #Calculate the propogation direction point
    
    Rotated_BeamProp=CoMTransposed_Rotated_PropPoint-CoMTransposed_Rotated_BeamFocusPoint 
    
    Rotated_BeamFocusPoint=CoMTransposed_Rotated_BeamFocusPoint+CentreOfMass
    
    return(Rotated_BeamFocusPoint,Rotated_BeamProp)
    
BeamFocus=np.array([[0],[0],[-2]])
BeamProp=np.array([[0],[0],[1]])
CentreOfMass=np.array([[0],[0],[1]])

New_BeamFocus,NewProp=RotateBeam(BeamFocus,BeamProp,CentreOfMass,0,-((np.pi)/2),0)

print(New_BeamFocus)
print(NewProp)
    