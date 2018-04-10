# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 21:55:35 2018

@author: Marius
"""

import numpy as np
import glob
import os
import subprocess
import time
import shutil
import random
import scipy.constants as constants

def ADDAForceConversion(Force,ElectricFieldStrength):

	Force=Force*(((ElectricFieldStrength)**2)/((constants.c)**2))*(10**(-5))
	
	return(Force)
	
def OurForceConversion(Force,CorrectionFactor,ElectricFieldStrength):
    
   return(ADDAForceConversion((Force*CorrectionFactor),ElectricFieldStrength))
   
def ElectricFieldStrengthCalc(DielConstant,Power,BeamWidth): #Power in Watts, BeamWidth in micro m
    
    Impedence=((constants.mu_0)/((constants.epsilon_0)*(DielConstant)))**0.5
    ElectricFieldStrength=((Impedence*Power)/((constants.pi)*(((BeamWidth*(1e-6))/2)**2)))**0.5
    
    return(ElectricFieldStrength)
    
def ReverseForce(Force, EFS, CF):
    return (Force*(constants.c**2)*(1e5))/(CF*(EFS**2))
    
    
CorrectionFactor=0.083861067
BeamWidth=1 #In micro m
Temperature=20 #Degrees C
Power=5e-3 #In Watts
MediumDielectricConstant=87.740-(0.40008*Temperature)+(9.398e-4*(Temperature**2))-(1.410e-6*(Temperature**3))
ElectricFieldStrength=ElectricFieldStrengthCalc(MediumDielectricConstant,Power,BeamWidth) #V/m

NewPower=2.5e-3

ForceFile = np.loadtxt("TotalForce")
OldForce = ReverseForce(ForceFile[:,1], ElectricFieldStrength, CorrectionFactor)
NEFS = ElectricFieldStrengthCalc(MediumDielectricConstant, NewPower, BeamWidth) #V/m
NewForce = OurForceConversion(OldForce, CorrectionFactor, NEFS)
ForceFile[:,1] = NewForce
np.savetxt("TotalForcePowerHalf", ForceFile, fmt="%e", delimiter=' ')