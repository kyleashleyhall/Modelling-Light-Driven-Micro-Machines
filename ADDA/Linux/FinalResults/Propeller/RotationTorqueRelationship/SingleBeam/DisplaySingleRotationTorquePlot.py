# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d

TitleString=r"Torque at single beam focus for a rotating propellor"

ax = plt.gca()
Plot = np.transpose(np.loadtxt('Torques'))
HalfArraySize=(np.size(Plot,axis=1)//2)+1
Value_A=((np.amax(Plot[1,:HalfArraySize])+np.amin(Plot[1,:HalfArraySize]))/2)
Value_B=((np.amax(Plot[1,:HalfArraySize])-np.amin(Plot[1,:HalfArraySize]))/2)
Value_Beta=90-(2*(Plot[0,(np.argmax(Plot[1,:HalfArraySize]))]))
Value_Sine=np.deg2rad(((2*Plot[0,:])+Value_Beta))
ForcePlotLabel=r"Raw Data; $A$="+str(Value_A)+r"; $B$="+str(Value_B)+r"; $\beta$="+str(Value_Beta)
ForcePlot=plt.errorbar(Plot[0,:], Plot[1,:],label=ForcePlotLabel)
PredictionPlotLabel=r"$A+B\sin(2\alpha + \beta)$"
PredictionPlot=plt.errorbar(Plot[0,:],(Value_A+(Value_B*(np.sin(Value_Sine)))),label=PredictionPlotLabel)
plt.title(TitleString)
plt.xlabel(r"Rotation angle, $\alpha$ (degrees)")
plt.ylabel(r"Torque (Nm)")
plt.legend(handles=[ForcePlot,PredictionPlot])
plt.show()