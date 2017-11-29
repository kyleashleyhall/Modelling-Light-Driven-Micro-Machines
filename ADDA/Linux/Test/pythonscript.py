# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



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
    
# Import Dipole Polorisations    
DipPolar = np.loadtxt('DipPol-Y', skiprows=1) #skip first row if it is a heading
x, y, z = DipPolar[:,0], DipPolar[:,1], DipPolar[:,2]
DipPolx = DipPolar[:,4] + 1j*DipPolar[:,5]
DipPoly = DipPolar[:,6] + 1j*DipPolar[:,7]
DipPolz = DipPolar[:,8] + 1j*DipPolar[:,9]

# Import Internal Electric Field contributions
Eintn = np.loadtxt('IntField-Y', skiprows=1)
Eintx = Eintn[:,4] + 1j*Eintn[:,5]
Einty = Eintn[:,6] + 1j*Eintn[:,7]
Eintz = Eintn[:,8] + 1j*Eintn[:,9]

# Import Electric field due to beam
Ebeam = np.loadtxt('IncBeam-Y', skiprows=1)
Ebeax = Ebeam[:,4] + 1j*Ebeam[:,5]
Ebeay = Ebeam[:,6] + 1j*Ebeam[:,7]
Ebeaz = Ebeam[:,8] + 1j*Ebeam[:,9]

# Import Dipole Seperation
DipSep = (np.loadtxt('DipSep'))+0 # To remove the array condition we add zero

axis = np.vstack([x,y,z]) 
DipolePolarisations = np.vstack([DipPolx, DipPoly, DipPolz]) #Dipole Polarisations
ElectricFields = np.vstack([(Eintx+Ebeax), (Einty+Ebeay), (Eintz+Ebeaz)]) #Internal Fields

CalculatedForces = Forces(axis, DipolePolarisations, ElectricFields, DipSep)
np.savetxt('Forces',np.transpose([np.vstack([axis,CalculatedForces])]),delimiter=' ')

plt.quiver(axis[0], axis[1], CalculatedForces[0], CalculatedForces[1])
plt.show()


print(sum(CalculatedForces[0]), sum(CalculatedForces[1]), sum(CalculatedForces[2]))











