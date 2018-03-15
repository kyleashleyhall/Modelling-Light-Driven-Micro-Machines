# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.mplot3d import axes3d

TitleString=r"Torque at single beam focus"

ax = plt.gca()
SpherePlot = np.transpose(np.loadtxt('..\..\..\Sphere\RotationTorqueRelationship\SingleBeam\Torques'))
SphereForcePlot=plt.errorbar(SpherePlot[0], SpherePlot[1],label="Sphere")
ConePlot = np.transpose(np.loadtxt('..\..\..\Cone\RotationTorqueRelationship\SingleBeam\Torques'))
ConeForcePlot=plt.errorbar(ConePlot[0], ConePlot[1],label="Cone")
PropellerPlot = np.transpose(np.loadtxt('..\..\..\Propeller\RotationTorqueRelationship\SingleBeam\Torques'))
PropellerForcePlot=plt.errorbar(PropellerPlot[0], PropellerPlot[1],label="Propeller")
plt.title(TitleString)
plt.xlabel(r"Rotation angle (degrees)")
plt.ylabel(r"Torque (Nm)")
plt.legend(handles=[SphereForcePlot,ConeForcePlot,PropellerForcePlot])
plt.show()