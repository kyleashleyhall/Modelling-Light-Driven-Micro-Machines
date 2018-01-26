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

def FileSlice(Fname,Ftype):
    filename = np.loadtxt(Fname, skiprows=1)	
    x,y,z = filename[:,0], filename[:,1], filename[:,2]
    if (Ftype=='float'):
        varX = filename[:,4]
        varY = filename[:,5]
        varZ = filename[:,6]
    elif (Ftype=='complex'):    
        varX = filename[:,4] + 1j*filename[:,5]
        varY = filename[:,6] + 1j*filename[:,7]
        varZ = filename[:,8] + 1j*filename[:,9]
    return np.vstack([x,y,z]), np.vstack([varX,varY,varZ])

def Grid(PositionsArray,GridType):
    Sep = DipSep(PositionsArray[0,:]) #dipsep only works in x axis
    GridSize = np.zeros([3], dtype=int) #Dimensions of Grid in x,y,z
    for i in range(3):
        minLocation = min(PositionsArray[i,:])
        maxLocation = max(PositionsArray[i,:])
        GridSize[i] = np.round((maxLocation-minLocation)/Sep)+1
    if (GridType=='float'):
        grid = np.zeros([GridSize[0],GridSize[1],GridSize[2],3], dtype=float)
        grid[:] = np.nan
    elif (GridType=='complex'):
        grid = np.zeros([GridSize[0],GridSize[1],GridSize[2],3], dtype=complex)
        grid[:] = np.nan
        
    return grid

def Forces(Axes, EField, DipoleSep, DipPol):
    
    Force=EField
    TotalForce=np.zeros([3])
    for i in range(len(EField[0])): #Iterate through all x values
        for j in range(len(EField[1])): #Iterate through all y values
            for k in range(len(EField[2])): #Iterate through all z values
                if (not np.isnan(Axes[i,j,k,0])): #Checks to see if the value is a dipole
                    
                    #Calculating the force (F_x)
                
                    CentralDifference=np.ones([2],dtype=bool) #To check if lower and upper part of Central difference is available
                    Force[i,j,k,0] = 0 #Initialise the value
                    try: #Check if there are x values below
                        if (np.isnan(Axes[(i-1),j,k,0])):
                            CentralDifference[0]=False
                    except:
                        CentralDifference[0]=False #If there are no x values below
                    try: #Check if there are x values above
                        if (np.isnan(Axes[(i+1),j,k,0])):
                            CentralDifference[1]=False
                    except:
                        CentralDifference[1]=False #If there are no x values above
                    if ((CentralDifference[0]==True)&(CentralDifference[1]==True)):
                        Force[i,j,k,0] += DipPol[i,j,k,0]*(np.conjugate((EField[(i+1),j,k,0]-EField[(i-1),j,k,0])/(2*DipoleSep)))
                        Force[i,j,k,0] += DipPol[i,j,k,1]*(np.conjugate((EField[(i+1),j,k,1]-EField[(i-1),j,k,1])/(2*DipoleSep)))
                        Force[i,j,k,0] += DipPol[i,j,k,2]*(np.conjugate((EField[(i+1),j,k,2]-EField[(i-1),j,k,2])/(2*DipoleSep)))
                        Force[i,j,k,0] = 0.5*np.real(Force[i,j,k,0])
                    elif ((CentralDifference[0]==False)&(CentralDifference[1]==True)):
                        Force[i,j,k,0] += DipPol[i,j,k,0]*(np.conjugate((EField[(i+1),j,k,0]-EField[i,j,k,0])/(DipoleSep)))
                        Force[i,j,k,0] += DipPol[i,j,k,1]*(np.conjugate((EField[(i+1),j,k,1]-EField[i,j,k,1])/(DipoleSep)))
                        Force[i,j,k,0] += DipPol[i,j,k,2]*(np.conjugate((EField[(i+1),j,k,2]-EField[i,j,k,2])/(DipoleSep)))
                        Force[i,j,k,0] = 0.5*np.real(Force[i,j,k,0])
                    elif ((CentralDifference[0]==True)&(CentralDifference[1]==False)):
                        Force[i,j,k,0] += DipPol[i,j,k,0]*(np.conjugate((EField[i,j,k,0]-EField[(i-1),j,k,0])/(DipoleSep)))
                        Force[i,j,k,0] += DipPol[i,j,k,1]*(np.conjugate((EField[i,j,k,1]-EField[(i-1),j,k,1])/(DipoleSep)))
                        Force[i,j,k,0] += DipPol[i,j,k,2]*(np.conjugate((EField[i,j,k,2]-EField[(i-1),j,k,2])/(DipoleSep)))
                        Force[i,j,k,0] = 0.5*np.real(Force[i,j,k,0])
                    else:
                        Force[i,j,k,0] = 0
        
                    #Calculating the force (F_y)

                    CentralDifference=np.ones([2],dtype=bool) #To check if lower and upper part of Central difference is available
                    Force[i,j,k,1] = 0 #Initialise the value
                    try: #Check if there are y values below
                        if (np.isnan(EField[i,(j-1),k,0])):
                            CentralDifference[0]=False
                    except:
                        CentralDifference[0]=False #If there are no y values below
                    try: #Check if there are y values above
                        if (np.isnan(EField[i,(j+1),k,0])):
                            CentralDifference[1]=False
                    except:
                        CentralDifference[1]=False #If there are no y values above
                    if ((CentralDifference[0]==True)&(CentralDifference[1]==True)):
                        Force[i,j,k,1] += DipPol[i,j,k,0]*(np.conjugate((EField[i,(j+1),k,0]-EField[i,(j-1),k,0])/(2*DipoleSep)))
                        Force[i,j,k,1] += DipPol[i,j,k,1]*(np.conjugate((EField[i,(j+1),k,1]-EField[i,(j-1),k,1])/(2*DipoleSep)))
                        Force[i,j,k,1] += DipPol[i,j,k,2]*(np.conjugate((EField[i,(j+1),k,2]-EField[i,(j-1),k,2])/(2*DipoleSep)))
                        Force[i,j,k,1] = 0.5*np.real(Force[i,j,k,1])
                    elif ((CentralDifference[0]==False)&(CentralDifference[1]==True)):
                        Force[i,j,k,1] += DipPol[i,j,k,0]*(np.conjugate((EField[i,(j+1),k,0]-EField[i,j,k,0])/(DipoleSep)))
                        Force[i,j,k,1] += DipPol[i,j,k,1]*(np.conjugate((EField[i,(j+1),k,1]-EField[i,j,k,1])/(DipoleSep)))
                        Force[i,j,k,1] += DipPol[i,j,k,2]*(np.conjugate((EField[i,(j+1),k,2]-EField[i,j,k,2])/(DipoleSep)))
                        Force[i,j,k,1] = 0.5*np.real(Force[i,j,k,1])
                    elif ((CentralDifference[0]==True)&(CentralDifference[1]==False)):
                        Force[i,j,k,1] += DipPol[i,j,k,0]*(np.conjugate((EField[i,j,k,0]-EField[i,(j-1),k,0])/(DipoleSep)))
                        Force[i,j,k,1] += DipPol[i,j,k,1]*(np.conjugate((EField[i,j,k,1]-EField[i,(j-1),k,1])/(DipoleSep)))
                        Force[i,j,k,1] += DipPol[i,j,k,2]*(np.conjugate((EField[i,j,k,2]-EField[i,(j-1),k,2])/(DipoleSep)))
                        Force[i,j,k,1] = 0.5*np.real(Force[i,j,k,1])
                    else:
                        Force[i,j,k,1] = 0
    
                    #Calculating the force (F_z)

                    CentralDifference=np.ones([2],dtype=bool) #To check if lower and upper part of Central difference is available
                    Force[i,j,k,2] = 0 #Initialise the value
                    try: #Check if there are z values below
                        if (np.isnan(EField[i,j,(k-1),0])):
                            CentralDifference[0]=False
                    except:
                        CentralDifference[0]=False #If there are no z values below
                    try: #Check if there are z values above
                        if (np.isnan(EField[i,j,(k+1),0])):
                            CentralDifference[1]=False
                    except:
                        CentralDifference[1]=False #If there are no z values above
                    if ((CentralDifference[0]==True)&(CentralDifference[1]==True)):
                        Force[i,j,k,2] += DipPol[i,j,k,0]*(np.conjugate((EField[i,j,(k+1),0]-EField[i,j,(k-1),0])/(2*DipoleSep)))
                        Force[i,j,k,2] += DipPol[i,j,k,1]*(np.conjugate((EField[i,j,(k+1),1]-EField[i,j,(k-1),1])/(2*DipoleSep)))
                        Force[i,j,k,2] += DipPol[i,j,k,2]*(np.conjugate((EField[i,j,(k+1),2]-EField[i,j,(k-1),2])/(2*DipoleSep)))
                        Force[i,j,k,2] = 0.5*np.real(Force[i,j,k,2])
                    elif ((CentralDifference[0]==False)&(CentralDifference[1]==True)):
                        Force[i,j,k,2] += DipPol[i,j,k,0]*(np.conjugate((EField[i,j,(k+1),0]-EField[i,j,k,0])/(DipoleSep)))
                        Force[i,j,k,2] += DipPol[i,j,k,1]*(np.conjugate((EField[i,j,(k+1),1]-EField[i,j,k,1])/(DipoleSep)))
                        Force[i,j,k,2] += DipPol[i,j,k,2]*(np.conjugate((EField[i,j,(k+1),2]-EField[i,j,k,2])/(DipoleSep)))
                        Force[i,j,k,2] = 0.5*np.real(Force[i,j,k,2])
                    elif ((CentralDifference[0]==True)&(CentralDifference[1]==False)):
                        Force[i,j,k,2] += DipPol[i,j,k,0]*(np.conjugate((EField[i,j,k,0]-EField[i,j,(k-1),0])/(DipoleSep)))
                        Force[i,j,k,2] += DipPol[i,j,k,1]*(np.conjugate((EField[i,j,k,1]-EField[i,j,(k-1),1])/(DipoleSep)))
                        Force[i,j,k,2] += DipPol[i,j,k,2]*(np.conjugate((EField[i,j,k,2]-EField[i,j,(k-1),2])/(DipoleSep)))
                        Force[i,j,k,2] = 0.5*np.real(Force[i,j,k,2])
                    else:
                        Force[i,j,k,2] = 0
                        
                    TotalForce[0]=TotalForce[0]+Force[i,j,k,0]
                    TotalForce[1]=TotalForce[1]+Force[i,j,k,1]
                    TotalForce[2]=TotalForce[2]+Force[i,j,k,2]
                    
                else:
                    Force[i,j,k,0]=np.nan
                    Force[i,j,k,1]=np.nan
                    Force[i,j,k,2]=np.nan
    
    print(TotalForce[0])
    print(TotalForce[1])
    print(TotalForce[2])
    
    return Force

#Preliminary Variables
Initial_dpl=30
Final_dpl=31
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
    DipPolRawPositions = FileSlice(DipFiles,'complex')[0] #Array containing the position of dipoles in DipPol
    DipPolRawValues = FileSlice(DipFiles,'complex')[1] #array with the polarization of said dipoles
    DipPolGridValues = Grid(DipPolRawPositions,'complex')
    AxesGrid = Grid(DipPolRawPositions,'float')
    IntFieldRawPositions = FileSlice(IntFFiles,'complex')[0] #Array containing the position of dipoles in IntField
    IntFieldRawValues = FileSlice(IntFFiles,'complex')[1] #array with the Internal field at each point
    IntFieldGridValues = Grid(IntFieldRawPositions,'complex')
    BeamRawPositions = FileSlice(BeamFiles,'complex')[0] #Array containing the position of dipoles in IncBeam
    BeamRawValues = FileSlice(BeamFiles,'complex')[1] #array with the Beam field at each point
    BeamGridValues = Grid(BeamRawPositions,'complex')
    RadForceRawPositions = FileSlice(ForceFiles,'float')[0] #Array containing the position of dipoles in RadForce
    RadForceRawValues = FileSlice(ForceFiles,'float')[1] #array with the ADDA Calculated RadForce at each point
    RadForceGridValues = Grid(RadForceRawPositions,'float')
    DipSeperation=DipSep(DipPolRawPositions[0,:]) #The Dipole seperation
    min_xLocation=min(DipPolRawPositions[0,:]) #Min x
    min_yLocation=min(DipPolRawPositions[1,:]) #Min y
    min_zLocation=min(DipPolRawPositions[2,:]) #Min z
    
    #Assign the DipPol Raw values to the DipPol Grid (including the axes positions)
    for i in range(len(DipPolRawPositions[1])):
        xLocation=np.int(np.round((DipPolRawPositions[0,i]-min_xLocation)/DipSeperation)) #Find the x Location in the grid
        yLocation=np.int(np.round((DipPolRawPositions[1,i]-min_yLocation)/DipSeperation)) #Find the y Location in the grid
        zLocation=np.int(np.round((DipPolRawPositions[2,i]-min_zLocation)/DipSeperation)) #Find the z Location in the grid
        AxesGrid[xLocation,yLocation,zLocation,0]=DipPolRawPositions[0,i]
        AxesGrid[xLocation,yLocation,zLocation,1]=DipPolRawPositions[1,i]
        AxesGrid[xLocation,yLocation,zLocation,2]=DipPolRawPositions[2,i]
        DipPolGridValues[xLocation,yLocation,zLocation,0]=DipPolRawValues[0,i]
        DipPolGridValues[xLocation,yLocation,zLocation,1]=DipPolRawValues[1,i]
        DipPolGridValues[xLocation,yLocation,zLocation,2]=DipPolRawValues[2,i]
        
    #Assign the IntField Raw values to the IntField Grid
    for i in range(len(IntFieldRawPositions[1])):
        xLocation=np.int(np.round((IntFieldRawPositions[0,i]-min_xLocation)/DipSeperation)) #Find the x Location in the grid
        yLocation=np.int(np.round((IntFieldRawPositions[1,i]-min_yLocation)/DipSeperation)) #Find the y Location in the grid
        zLocation=np.int(np.round((IntFieldRawPositions[2,i]-min_zLocation)/DipSeperation)) #Find the z Location in the grid
        IntFieldGridValues[xLocation,yLocation,zLocation,0]=IntFieldRawValues[0,i]
        IntFieldGridValues[xLocation,yLocation,zLocation,1]=IntFieldRawValues[1,i]
        IntFieldGridValues[xLocation,yLocation,zLocation,2]=IntFieldRawValues[2,i]
        
    #Assign the IncBeam Raw value to the IncBeam Grid
    for i in range(len(BeamRawPositions[1])):
        xLocation=np.int(np.round((BeamRawPositions[0,i]-min_xLocation)/DipSeperation)) #Find the x Location in the grid
        yLocation=np.int(np.round((BeamRawPositions[1,i]-min_yLocation)/DipSeperation)) #Find the y Location in the grid
        zLocation=np.int(np.round((BeamRawPositions[2,i]-min_zLocation)/DipSeperation)) #Find the z Location in the grid
        BeamGridValues[xLocation,yLocation,zLocation,0]=BeamRawValues[0,i]
        BeamGridValues[xLocation,yLocation,zLocation,1]=BeamRawValues[1,i]
        BeamGridValues[xLocation,yLocation,zLocation,2]=BeamRawValues[2,i]
        
    #Assign the RadForce values to the RadFroce Grid
    for i in range(len(RadForceRawPositions[1])):
        xLocation=np.int(np.round((RadForceRawPositions[0,i]-min_xLocation)/DipSeperation)) #Find the x Location in the grid
        yLocation=np.int(np.round((RadForceRawPositions[1,i]-min_yLocation)/DipSeperation)) #Find the y Location in the grid
        zLocation=np.int(np.round((RadForceRawPositions[2,i]-min_zLocation)/DipSeperation)) #Find the z Location in the grid
        RadForceGridValues[xLocation,yLocation,zLocation,0]=RadForceRawValues[0,i]
        RadForceGridValues[xLocation,yLocation,zLocation,1]=RadForceRawValues[1,i]
        RadForceGridValues[xLocation,yLocation,zLocation,2]=RadForceRawValues[2,i]
    
    EFieldGridValues = IntFieldGridValues + BeamGridValues
    CalcForce = Forces(AxesGrid,EFieldGridValues,DipSeperation,DipPolGridValues)
    EndTime=time.clock()
    
    #SAVE CALCULATED FORCES
    with open(FFiles,'wb') as f:
        PrintedForce=np.zeros([len(RadForceRawPositions[1]),7])
        Iterator=0
        f.write(b'x y z |F|^2 Fx Fy Fz \n')
        EstimatedParticleForce=np.zeros([3])
        for i in range(len(CalcForce[0])): #Iterate through all x values
            for j in range(len(CalcForce[1])): #Iterate through all y values
                for k in range(len(CalcForce[2])): #Iterate through all z values
                    if (not np.isnan(AxesGrid[i,j,k,0])): #Checks to see if the value is a dipole
                        PrintedForce[Iterator,0]=AxesGrid[i,j,k,0]
                        PrintedForce[Iterator,1]=AxesGrid[i,j,k,1]
                        PrintedForce[Iterator,2]=AxesGrid[i,j,k,2]
                        PrintedForce[Iterator,3]=pow((pow(np.real(CalcForce[i,j,k,0]),2)+pow(np.real(CalcForce[i,j,k,1]),2)+pow(np.real(CalcForce[i,j,k,2]),2)),0.5)
                        PrintedForce[Iterator,4]=np.real(CalcForce[i,j,k,0])
                        PrintedForce[Iterator,5]=np.real(CalcForce[i,j,k,1])
                        PrintedForce[Iterator,6]=np.real(CalcForce[i,j,k,2])
                        Iterator += 1
        
        np.savetxt(f,PrintedForce, fmt='%.10f',delimiter=' ')
    
    #This section is where we look at the ADDA Calculated Forces
    EstimatedParticleForce=np.array([[np.sum(PrintedForce[:,4])],[np.sum(PrintedForce[:,5])],[np.sum(PrintedForce[:,6])]])
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
    
    #Use to delete the files after processing
    '''try:
        shutil.rmtree(FFiles.replace(os.sep+'CalculatedForces',''))
    except:
        print('Cannot Delete')'''
    dpl=dpl+Step_dpl

ForceErrorPath = str(os.getcwd())+str(os.sep+'ForceError') #|(F_ADDA - F_Calc)|/F_ADDA	
with open(ForceErrorPath, 'wb') as f:
    f.write(b'dpl |F_ADDA(x)-F_Calc(x)| |F_ADDA(y)-F_Calc(y)| |F_ADDA(z)-F_Calc(z)| |F_ADDA(x)-F_Calc(x)|/F_ADDA(x) |F_ADDA(y)-F_Calc(y)|/F_ADDA(y) |F_ADDA(z)-F_Calc(z)|/F_ADDA(z)\n')
    np.savetxt(f, ForceError, fmt='%.10f', delimiter=' ')
TimeLogPath = str(os.getcwd())+str(os.sep+'TimeLog')	
with open(TimeLogPath, 'wb') as f:
    f.write(b'dpl Time(s)\n')
    np.savetxt(f, TimeRecordings, fmt='%.10f', delimiter=' ')