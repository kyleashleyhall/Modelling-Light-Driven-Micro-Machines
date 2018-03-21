# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d

TitleString=r"Cone torque at single beam focus"

ax = plt.gca()
'''File1 = np.transpose(np.loadtxt('DipoleSep_0_03\Torques'))
Plot1=plt.errorbar(File1[0], File1[1],label=r"Dipole Separation 0.03$\mu m$")'''
File2 = np.transpose(np.loadtxt('DipoleSep_0_04\Torques'))
Plot2=plt.errorbar(File2[0], File2[1],label=r"Dipole Separation 0.04$\mu m$")
File3 = np.transpose(np.loadtxt('DipoleSep_0_05\Torques'))
Plot3=plt.errorbar(File3[0], File3[1],label=r"Dipole Separation 0.05$\mu m$")
File4 = np.transpose(np.loadtxt('DipoleSep_0_06\Torques'))
Plot4=plt.errorbar(File4[0], File4[1],label=r"Dipole Separation 0.06$\mu m$")
File5 = np.transpose(np.loadtxt('DipoleSep_0_07\Torques'))
Plot5=plt.errorbar(File5[0], File5[1],label=r"Dipole Separation 0.07$\mu m$")
File6 = np.transpose(np.loadtxt('DipoleSep_0_08\Torques'))
Plot6=plt.errorbar(File6[0], File6[1],label=r"Dipole Separation 0.08$\mu m$")
File7 = np.transpose(np.loadtxt('DipoleSep_0_09\Torques'))
Plot7=plt.errorbar(File7[0], File7[1],label=r"Dipole Separation 0.09$\mu m$")
plt.title(TitleString)
plt.xlabel(r"Rotation angle (degrees)")
plt.ylabel(r"Torque (Nm)")
#plt.legend(handles=[Plot1,Plot2,Plot3,Plot4,Plot5,Plot6,Plot7])
plt.legend(handles=[Plot2,Plot3,Plot4,Plot5,Plot6,Plot7])
plt.show()