# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants

def ADDAForceConversion(Force,ElectricFieldStrength):

	Force=Force*(((ElectricFieldStrength)**2)/((constants.c)**2))*(10**(-5))
	
	return(Force)
    
def OurForceConversion(Force,CorrectionFactor,ElectricFieldStrength):
    
   return(ADDAForceConversion((Force*CorrectionFactor),ElectricFieldStrength))
   
ElectricFieldStrength=1 #V/m
CorrectionFactor=0.308310573308


Forces = np.transpose(np.loadtxt('TotalForce'))

Forces[3]=OurForceConversion(Forces[3],CorrectionFactor,ElectricFieldStrength)
Forces[4]=OurForceConversion(Forces[4],CorrectionFactor,ElectricFieldStrength)
Forces[5]=OurForceConversion(Forces[5],CorrectionFactor,ElectricFieldStrength)

plt.figure(1)

plt.subplot(311)
plt.plot(Forces[2],Forces[3])

plt.subplot(312)
plt.plot(Forces[2],Forces[4])

plt.subplot(313)
plt.plot(Forces[2],Forces[5])

plt.show()