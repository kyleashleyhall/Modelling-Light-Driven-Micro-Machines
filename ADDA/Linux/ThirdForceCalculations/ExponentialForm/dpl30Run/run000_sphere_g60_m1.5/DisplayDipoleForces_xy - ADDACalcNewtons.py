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


ElectricFieldStrength=1 #V/m

Quiv = np.transpose(np.loadtxt('RadForce-Y', skiprows=1))

Quiv[4]=ADDAForceConversion(Quiv[4],ElectricFieldStrength)
Quiv[5]=ADDAForceConversion(Quiv[5],ElectricFieldStrength)
Quiv[6]=ADDAForceConversion(Quiv[6],ElectricFieldStrength)

print(sum(Quiv[4]))
print(sum(Quiv[5]))
print(sum(Quiv[6]))

plt.quiver(Quiv[0], Quiv[1], Quiv[4], Quiv[5])
plt.show()