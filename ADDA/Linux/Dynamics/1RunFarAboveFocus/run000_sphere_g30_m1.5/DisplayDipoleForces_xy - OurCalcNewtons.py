# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import scipy.constants as constants

def ADDAForceConversion(Force,ElectricFieldStrength):

	Force=Force*(((ElectricFieldStrength)**2)/((constants.c)**2))*(10**(-5))
	
	return(Force)
    
def OurForceConversion(Force,CorrectionFactor,ElectricFieldStrength):
    
   return(ADDAForceConversion((Force*CorrectionFactor),ElectricFieldStrength))
   
ElectricFieldStrength=1 #V/m
CorrectionFactor=0.308310573308

Quiv = np.transpose(np.loadtxt('CalculatedForces', skiprows=1))

Quiv[4]=OurForceConversion(Quiv[4],CorrectionFactor,ElectricFieldStrength)
Quiv[5]=OurForceConversion(Quiv[5],CorrectionFactor,ElectricFieldStrength)
Quiv[6]=OurForceConversion(Quiv[6],CorrectionFactor,ElectricFieldStrength)

print(sum(Quiv[4]))
print(sum(Quiv[5]))
print(sum(Quiv[6]))

plt.quiver(Quiv[0], Quiv[1], Quiv[4], Quiv[5])
plt.show()