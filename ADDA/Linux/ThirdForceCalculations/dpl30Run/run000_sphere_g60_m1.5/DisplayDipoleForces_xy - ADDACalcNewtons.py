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


Quiv = np.transpose(np.loadtxt('RadForce-Y', skiprows=1))

print(ADDAForceConversion(sum(Quiv[4]),1))
print(ADDAForceConversion(sum(Quiv[5]),1))
print(ADDAForceConversion(sum(Quiv[6]),1))

plt.quiver(Quiv[0], Quiv[1], ADDAForceConversion(Quiv[4],1), ADDAForceConversion(Quiv[5],1))
plt.show()