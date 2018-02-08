# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius and Kyle
"""

import numpy as np

ADDAForces = np.transpose(np.loadtxt('RadForce-Y', skiprows=1))

#Either skiprows=1 or skiprows=0 depending on the maturity of the project
OurForces = np.transpose(np.loadtxt('CalculatedForces', skiprows=1))

print((sum(ADDAForces[6]))/(sum(OurForces[6])))
input()