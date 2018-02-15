# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:50:50 2018

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

def CalculateCentreOfMass(xPositions,yPositions,zPositions):
    
    x_CoM=(np.sum(xPositions))/(np.size(xPositions))
    y_CoM=(np.sum(yPositions))/(np.size(yPositions))
    z_CoM=(np.sum(zPositions))/(np.size(zPositions))
    
    CentreOfMass=np.array([[x_CoM],[y_CoM],[z_CoM]])
    
    return(CentreOfMass)
    
def CalculateMomentOfInertia(xPositions,yPositions,zPositions,CentreOfMass,Mass):
    
    DipMass=Mass/(np.size(xPositions))
    Inertia=np.zeros([3,1])
    Inertia[0,0]=np.sum((DipMass*((np.square(((yPositions-CentreOfMass[1,0])*(10**(-6)))))+(np.square(((zPositions-CentreOfMass[2,0])*(10**(-6))))))))
    Inertia[1,0]=np.sum((DipMass*((np.square(((xPositions-CentreOfMass[0,0])*(10**(-6)))))+(np.square(((zPositions-CentreOfMass[2,0])*(10**(-6))))))))
    Inertia[2,0]=np.sum((DipMass*((np.square(((xPositions-CentreOfMass[0,0])*(10**(-6)))))+(np.square(((yPositions-CentreOfMass[1,0])*(10**(-6))))))))
    
    return(Inertia)
    
def CalculateTorque(xPositions,yPositions,zPositions,F_x,F_y,F_z): #Assumes Forces are in Netwons!

    Torque=np.zeros([6,np.size(xPositions)]) #Create an empty array to store position and the Torque in Nm at that point
    
    for i in range(np.size(xPositions)):
        Torque[0,i]=xPositions[i]
        Torque[1,i]=yPositions[i]
        Torque[2,i]=zPositions[i]
        Torque[3,i]=(((yPositions[i])*(F_z[i]))-((zPositions[i])*(F_y[i])))
        Torque[4,i]=(((zPositions[i])*(F_x[i]))-((xPositions[i])*(F_z[i])))
        Torque[5,i]=(((xPositions[i])*(F_y[i]))-((yPositions[i])*(F_x[i])))
        
    return(Torque)

def CalculateAngles(Omega_0,TotalTorque,Inertia,TimeStep):
    
    AngleDisplacement=(Omega_0*TimeStep)+(0.5*(TotalTorque/Inertia)*(TimeStep**2)) #In Rad
    NewOmega=Omega_0+((TotalTorque/Inertia)*TimeStep) #In Rads^-1
    
    return(AngleDisplacement,NewOmega)

CalculatedForces=np.transpose(np.loadtxt('CalculatedForces', skiprows=1))
Mass=1
TimeStep=5e-5

CentreOfMass=CalculateCentreOfMass(CalculatedForces[0,:],CalculatedForces[1,:],CalculatedForces[2,:])
Inertia=CalculateMomentOfInertia(CalculatedForces[0,:],CalculatedForces[1,:],CalculatedForces[2,:],CentreOfMass,Mass)
DipoleTorque=CalculateTorque(CalculatedForces[0,:],CalculatedForces[1,:],CalculatedForces[2,:],CalculatedForces[4,:],CalculatedForces[5,:],CalculatedForces[6,:])
Omega_theta,Omega_psi,Omega_phi=0,0,0
Theta,Omega_theta=CalculateAngles(Omega_theta,np.sum(DipoleTorque[3,:]),Inertia[0,0],TimeStep)
Psi,Omega_psi=CalculateAngles(Omega_psi,np.sum(DipoleTorque[4,:]),Inertia[1,0],TimeStep)
Phi,Omega_phi=CalculateAngles(Omega_phi,np.sum(DipoleTorque[5,:]),Inertia[2,0],TimeStep)

print(Theta)
print(Psi)
print(Phi)

print(Omega_theta)
print(Omega_psi)
print(Omega_phi)