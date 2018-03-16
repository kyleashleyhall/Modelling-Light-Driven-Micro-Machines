# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d

TitleString=r"Torque at single beam focus for a rotating cone"

ax = plt.gca()
Plot = np.transpose(np.loadtxt('Torques'))
ForcePlot=plt.plot(Plot[0], Plot[1])
plt.title(TitleString)
plt.xlabel(r"Rotation angle (degrees)")
plt.ylabel(r"Torque (Nm)")
plt.show()