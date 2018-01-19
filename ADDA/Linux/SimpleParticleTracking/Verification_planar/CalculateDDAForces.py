# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import numpy as np
import glob
import os
import subprocess
import time
import shutil
import fnmatch

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
				Force[j,i]=Force[j,i]+((0.5*(np.real((DipPol[k,i])*(Eint[k,i])*(np.conjugate(dE[k,j,i]))))))
	
	return Force
	
def DipSep(Singleaxis):
	dx = np.zeros([len(Singleaxis)-1])
	for i in range(len(Singleaxis)-1):
		dx[i] = Singleaxis[i+1] - Singleaxis[i]
	return min(np.absolute(dx))

#Preliminary Variables
Initial_dpl=15
Final_dpl=16
Step_dpl=1


#Perform the DDA Calculations and calculate forces
DipPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'DipPol-Y')  
IntFPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IntField-Y')
BeamPathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'IncBeam-Y')
ForcePathInput = str(os.getcwd())+str(os.sep+'*'+os.sep+'RadForce-Y')
dpl=Initial_dpl
ForceError=np.zeros([(Final_dpl-Initial_dpl),7])
TimeRecordings=np.zeros([(Final_dpl-Initial_dpl),2])
CalculationTimes=np.zeros([(Final_dpl-Initial_dpl),2])
while (dpl<Final_dpl):
	StartTime=time.clock()
	print('Processing dpl: '+str(dpl))
	callString=".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl "+str(dpl)+" -sym enf -lambda 1 -prop 0 0 1 -store_beam -store_dip_pol -store_int_field -store_force" #The script for performing the DDA calculations
	print(".."+os.sep+"src"+os.sep+"seq"+os.sep+"adda -size 2 -dpl "+str(dpl)+" -sym enf -lambda 1 -prop 0 0 1 -store_beam -store_dip_pol -store_int_field -store_force")
	subprocess.call(callString,shell=True)
	DipFiles, IntFFiles, BeamFiles, ForceFiles = sorted(glob.glob(DipPathInput))[-1], sorted(glob.glob(IntFPathInput))[-1], sorted(glob.glob(BeamPathInput))[-1], sorted(glob.glob(ForcePathInput))[-1] #File containing the paths to each DipPol, IntField file
	FFiles = DipFiles.replace('DipPol-Y','CalculatedForces')
	axes = FileSlice(DipFiles)[0] #array containing the position of all dipoles
	DipPol = FileSlice(DipFiles)[1] #array with the polarization of said dipoles
	IntField = FileSlice(IntFFiles)[1] #array with the Internal field at each point
	BeamF = FileSlice(BeamFiles)[1] #array with the Beam field at each point
	EField = IntField + BeamF #Total electric field 
	DipSeperation = DipSep(axes[0]) #Calculate the dipole seperation
	CalcForce = Forces(axes, DipPol, EField, DipSeperation)
	EndTime=time.clock()
	np.savetxt(FFiles,np.transpose([np.vstack([axes,CalcForce])]), fmt='%.10f',delimiter=' ')
	EstimatedParticleForce = np.array([[np.sum(CalcForce[0])],[np.sum(CalcForce[1])],[np.sum(CalcForce[2])]])
	
	#This section is where we look at the ADDA Calculated Forces
	ADDADipoleForceFile = np.loadtxt(ForceFiles, skiprows=1) #Load the ADDA Dipole Forces File
	ADDAParticleForce = np.array([[np.sum(ADDADipoleForceFile[:,4])],[np.sum(ADDADipoleForceFile[:,5])],[np.sum(ADDADipoleForceFile[:,6])]]) #Save the ADDA Particle Forces to memory
	ForceError[(dpl-Initial_dpl),0] = dpl 
	ForceError[(dpl-Initial_dpl),1] = pow((pow(ADDAParticleForce[0]-EstimatedParticleForce[0],2)),0.5)
	ForceError[(dpl-Initial_dpl),2] = pow((pow(ADDAParticleForce[1]-EstimatedParticleForce[1],2)),0.5)
	ForceError[(dpl-Initial_dpl),3] = pow((pow(ADDAParticleForce[2]-EstimatedParticleForce[2],2)),0.5)
	ForceError[(dpl-Initial_dpl),4] = pow((pow(ADDAParticleForce[0]-EstimatedParticleForce[0],2)),0.5)/ADDAParticleForce[0]
	ForceError[(dpl-Initial_dpl),5] = pow((pow(ADDAParticleForce[1]-EstimatedParticleForce[1],2)),0.5)/ADDAParticleForce[1]
	ForceError[(dpl-Initial_dpl),6] = pow((pow(ADDAParticleForce[2]-EstimatedParticleForce[2],2)),0.5)/ADDAParticleForce[2]
	TimeRecordings[(dpl-Initial_dpl),0] = dpl
	TimeRecordings[(dpl-Initial_dpl),1] = EndTime-StartTime
	try:
		shutil.rmtree(FFiles.replace(os.sep+'CalculatedForces',''))
	except:
		print('Cannot Delete')
	dpl=dpl+Step_dpl

ForceErrorPath = str(os.getcwd())+str(os.sep+'ForceError') #|(F_ADDA - F_Calc)|/F_ADDA	
with open(ForceErrorPath, 'wb') as f:
	f.write(b'dpl |F_ADDA(x)-F_Calc(x)| |F_ADDA(y)-F_Calc(y)| |F_ADDA(z)-F_Calc(z)| |F_ADDA(x)-F_Calc(x)|/F_ADDA(x) |F_ADDA(y)-F_Calc(y)|/F_ADDA(y) |F_ADDA(z)-F_Calc(z)|/F_ADDA(z)\n')
	#f.write(bytes("SP,"+lists+"\n","UTF-8"))
	#Used this line for a variable list of numbers
	np.savetxt(f, ForceError, fmt='%.10f', delimiter=' ')
TimeLogPath = str(os.getcwd())+str(os.sep+'TimeLog')	
with open(TimeLogPath, 'wb') as f:
	f.write(b'dpl Time(s)\n')
	#f.write(bytes("SP,"+lists+"\n","UTF-8"))
	#Used this line for a variable list of numbers
	np.savetxt(f, TimeRecordings, fmt='%.10f', delimiter=' ')










