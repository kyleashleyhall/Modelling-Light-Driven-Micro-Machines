# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 15:09:12 2018

@author: Marius
"""

import os
import numpy as np
import glob
import matplotlib.pyplot as plt

def potential(px, temp):
    bzc = 1.38064853e-23
    return - (temp*bzc) * np.log(px)

filepaths = glob.glob(str(os.getcwd())+os.sep+'halfpower'+os.sep+'*')
allfiles = np.zeros((0,3))
for i in range(len(filepaths)):
    file = np.loadtxt(filepaths[i])
    allfiles = np.vstack((allfiles, file))
    #plt.plot(file[:,0], file[:,1], label=str(filepaths[i][-1]))plt.show()

#plt.hist(allfiles[:,1], bins=100)
#plt.show()

binheight, binpositions = np.histogram(allfiles[:,1], bins=100, normed=True)
temp = 273+20
U = potential(binheight, temp)
plt.plot(binpositions[:-1], U)
plt.show()
