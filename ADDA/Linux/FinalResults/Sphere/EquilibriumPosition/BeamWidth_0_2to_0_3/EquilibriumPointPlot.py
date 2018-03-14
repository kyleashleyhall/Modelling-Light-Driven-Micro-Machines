# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 20:48:12 2018

@author: Marius
"""

import os
import numpy as np
import glob
import matplotlib.pyplot as plt

FilesLoc1 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_2to0_209'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc2 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_21to0_219'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc3 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_22to0_229'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc4 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_23to0_239'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc5 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_24to0_249'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc6 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_25to0_259'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc7 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_26to0_269'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc8 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_27to0_279'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc9 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_28to0_289'+os.sep+'*'+os.sep+'EquilibriumPosition'))
FilesLoc10 = (str(os.getcwd())+str(os.sep+'BeamWidth_0_29to0_3'+os.sep+'*'+os.sep+'EquilibriumPosition'))
EPositionsArray = np.zeros((0,2))
BeamWidth= np.linspace(0.2,0.3,101)
for i in range(10):
    BeamWidth= np.linspace(0.2,0.209,10)
    EqFile = (glob.glob(FilesLoc1))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)
    
for i in range(10):
    BeamWidth= np.linspace(0.21,0.219,10)
    EqFile = (glob.glob(FilesLoc2))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)

for i in range(10):
    BeamWidth= np.linspace(0.22,0.229,10)
    EqFile = (glob.glob(FilesLoc3))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)
    
for i in range(10):
    BeamWidth= np.linspace(0.23,0.239,10)
    EqFile = (glob.glob(FilesLoc4))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)
    
for i in range(10):
    BeamWidth= np.linspace(0.24,0.249,10)
    EqFile = (glob.glob(FilesLoc5))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)
    
for i in range(10):
    BeamWidth= np.linspace(0.25,0.259,10)
    EqFile = (glob.glob(FilesLoc6))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)
    
for i in range(10):
    BeamWidth= np.linspace(0.26,0.269,10)
    EqFile = (glob.glob(FilesLoc7))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)
    
for i in range(10):
    BeamWidth= np.linspace(0.27,0.279,10)
    EqFile = (glob.glob(FilesLoc8))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)
    
for i in range(10):
    BeamWidth= np.linspace(0.28,0.289,10)
    EqFile = (glob.glob(FilesLoc9))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)

for i in range(11):
    BeamWidth= np.linspace(0.29,0.3,11)
    EqFile = (glob.glob(FilesLoc10))[i]
    EqFile1 = np.loadtxt(EqFile, skiprows=1)
    EquilibriumValue = EqFile1[-1,1]
    EPositionsArray = np.append(EPositionsArray, np.array([[BeamWidth[i], EquilibriumValue]]),axis=0)
    
plt.plot(EPositionsArray[:,0], EPositionsArray[:,1])
plt.xlabel('Beamwidth')
plt.ylabel('Equilibrium position in z')
plt.show()
    