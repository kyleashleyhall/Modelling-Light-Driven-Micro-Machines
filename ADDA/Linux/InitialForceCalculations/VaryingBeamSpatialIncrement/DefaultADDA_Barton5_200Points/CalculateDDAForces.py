# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import numpy as np
import glob
import os
import subprocess
import shutil

def FileSlice(Fname):
    filename = np.loadtxt(Fname, skiprows=1)    
    x,y,z = filename[:,0], filename[:,1], filename[:,2]
    varX = filename[:,4] + 1j*filename[:,5]
    varY = filename[:,6] + 1j*filename[:,7]
    varZ = filename[:,8] + 1j*filename[:,9]
    return np.vstack([x,y,z]), np.vstack([varX,varY,varZ])

def Forces(axes, DipPol, Eint, DipoleSep):
    d_l = np.zeros([2])  #d for central difference formula, for each axis
    dE = np.zeros([3, 3, len(Eint[0])], dtype=complex) #Derivative at each point, for each axis
    Force = np.zeros([3, len(Eint[0])]) # Force at each point
    CurrentPosition = np.zeros([3])
    for i in range(len(axes[0])): #Iterates through each point in the file

      #Varying through each axis
        for j in range(3):
            ValuesAvailable=np.array([True,True]) #We assume that there are X,Y or Z Values available before and after
            if (j==0): #Calculating dx values
                IndexValues=(np.where((axes[1,:]==axes[1,i])&(axes[2,:]==axes[2,i])))[0]#Finds the row numbers of all values with the same Y,Z positions (varying X)
            elif (j==1): #Calculating dy values
                IndexValues=(np.where((axes[0,:]==axes[0,i])&(axes[2,:]==axes[2,i])))[0]#Finds the row numbers of all values with the same X,Z positions (varying Y)
            elif (j==2): #Calculating dz values
                IndexValues=(np.where((axes[1,:]==axes[1,i])&(axes[0,:]==axes[0,i])))[0]#Finds the row numbers of all values with the same Y,X positions (varying Z)
            try: #Attempt to find a value before
                d_l[0]=((np.where(IndexValues<i))[0])[-1] #The preceeding Ex, Ey and Ez values
            except: #When no value before can be found
                ValuesAvailable[0]=False #We have not found a preceeding value
            try: #Attempt to find a value after
                d_l[1]=((np.where(IndexValues>i))[0])[0] #The following Ex, Ey and Ez values
            except:
                ValuesAvailable[1]=False #We have not found a following value
	#Quick Sanity check to ensure the difference is reasonable:
	#Come Back To This
	#Calculate dE/dj
            if ((ValuesAvailable[0]==True)&(ValuesAvailable[1]==True)):
	  #Calculate dE_l/dj by the central difference theorem 
                for k in range(3):
                    dE[k,j,i]=((Eint[k,(d_l[1])]-Eint[k,(d_l[0])])/(DipoleSep*2))
            elif ((ValuesAvailable[0]==True)&(ValuesAvailable[1]==False)):
            #Calculate dE_l/dj by using only the preceeding value
                for k in range(3):
                    dE[k,j,i]=((Eint[k,i]-Eint[k,(d_l[0])])/(DipoleSep))
            elif ((ValuesAvailable[0]==False)&(ValuesAvailable[1]==True)):
	  #Calculate dE_l/dj by using only the following value
                for k in range(3):
                    dE[k,j,i]=((Eint[k,(d_l[1])]-Eint[k,i])/(DipoleSep))
            else:
	      #Set (dE_x)/(dj), (dE_y)/(dj), (dE_z)/(dj) equal to zero
                for k in range(3):
                    dE[k,j,i]=0

      #Calculating the force
        for j in range(3): #Force in each direction	
            for k in range(3): #Summating the different contributions
                Force[j,i]=Force[j,i]+(0.5*(np.real((DipPol[k,i])*(np.conjugate(dE[k,j,i])))))
    
    return Force
    
def DipSep(Singleaxis):
    dx = np.zeros([len(Singleaxis)-1])
    for i in range(len(Singleaxis)-1):
        dx[i] = Singleaxis[i+1] - Singleaxis[i]
    return min(np.absolute(dx))

#Preliminary Variables
x=0
z=5
Initial_y=-5
Step_y=0.05
Final_y=5

y=Initial_y #Set the value of y to the initial value


#Perform the DDA Calculations and calculate forces
DipPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'DipPol-Y')  
IntFPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IntField-Y')
BeamPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IncBeam-Y')
while (y<Final_y):
  print('\nProcessing DDA for beam location x='+str(x)+', y='+str(y)+' z='+str(z)+'...\n')
  callString=".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -sym enf -lambda 1 -prop 0 0 1 -beam barton5 1 "+str(x)+" "+str(y)+" "+str(z)+" -store_beam -store_dip_pol -store_int_field" #The script for performing the DDA calculations
  print(".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -sym enf -lambda 1 -prop 0 0 1 -beam barton5 1 "+str(x)+" "+str(y)+" "+str(z)+" -store_beam -store_dip_pol -store_int_field")
  subprocess.call(callString,shell=True)
  DipFiles, IntFFiles, BeamFiles = sorted(glob.glob(DipPathInput))[-1], sorted(glob.glob(IntFPathInput))[-1], sorted(glob.glob(BeamPathInput))[-1] #File containing the paths to each DipPol, IntField file
  FFiles = DipFiles.replace('DipPol-Y','Forces')
  axes = FileSlice(DipFiles)[0] #array containing the position of all dipoles
  DipPol = FileSlice(DipFiles)[1] #array with the polarization of said dipoles
  IntField = FileSlice(IntFFiles)[1] #array with the Internal field at each point
  BeamF = FileSlice(BeamFiles)[1] #array with the Beam field at each point
  EField = IntField + BeamF #Total electric field 
  DipSeperation = DipSep(axes[0]) #Calculate the dipole seperation
  print('\nProcessing Forces for beam location x='+str(x)+', y='+str(y)+' z='+str(z)+'...\n')
  CalcForce = Forces(axes, DipPol, EField, DipSeperation)
  np.savetxt(FFiles,np.transpose([np.vstack([axes,CalcForce])]), fmt='%.10f',delimiter=' ')
  ParticleForce = np.array([[x],[y],[z],[sum(CalcForce[0])],[sum(CalcForce[1])],[sum(CalcForce[2])]])
  try:
    shutil.rmtree(FFiles.replace(os.sep+'Forces',''))
  except:
    print('Cannot Delete')
  try:
    AllForces = np.hstack([AllForces,ParticleForce])
  except:
    AllForces = ParticleForce
  y=y+Step_y
AllForcesPath = str(os.getcwd())+str(os.sep+'TotalForce')
np.savetxt(AllForcesPath,np.transpose(AllForces), fmt='%.10f',delimiter=' ')










