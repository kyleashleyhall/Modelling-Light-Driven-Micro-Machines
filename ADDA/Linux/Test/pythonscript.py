# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:57:26 2017

@author: Marius
"""

import numpy as np
import matplotlib.pyplot as plt



rfile = np.loadtxt('DipPol-Y', skiprows=1) #skip first row if it is a heading
print(rfile[1,1])

