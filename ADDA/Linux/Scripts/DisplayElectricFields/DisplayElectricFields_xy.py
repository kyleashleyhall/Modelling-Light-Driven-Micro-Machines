# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 11:05:22 2018

@author: Marius and Kyle
"""

import matplotlib.pyplot as plt
import numpy as np

def DipSep(Singleaxis):
    dx = np.zeros([len(Singleaxis)-1])
    for i in range(len(Singleaxis)-1):
        dx[i] = Singleaxis[i+1] - Singleaxis[i]
        if (dx[i]==0):
            dx[i]=np.inf
    return min(np.absolute(dx))

FileName=input('Filename (IntField): ')

IntFieldRaw = np.transpose(np.loadtxt(FileName,skiprows=1))

FileName=input('Filename (IncBeam): ')

IncBeamRaw = np.transpose(np.loadtxt(FileName,skiprows=1))

DipoleSeperation=DipSep(IntFieldRaw[0,:])

#Value in x to fix:
Fixed_z=-0.08333169654

Fixed_zLocations_IntField=np.argwhere(((Fixed_z-(DipoleSeperation/2))<IntFieldRaw[2,:])&(IntFieldRaw[2,:]<(Fixed_z+(DipoleSeperation/2))))

axis_array = np.zeros([(len(Fixed_zLocations_IntField)),3], dtype=float)
EField_array = np.zeros([len(Fixed_zLocations_IntField),3], dtype=complex)

for i in range(len(Fixed_zLocations_IntField)):
    IncBeamLineNumber=np.argwhere(((IntFieldRaw[0,(Fixed_zLocations_IntField[i])]-(DipoleSeperation/2))<IncBeamRaw[0,:])&(IncBeamRaw[0,:]<(IntFieldRaw[0,(Fixed_zLocations_IntField[i])]+(DipoleSeperation/2)))&((IntFieldRaw[1,(Fixed_zLocations_IntField[i])]-(DipoleSeperation/2))<IncBeamRaw[1,:])&(IncBeamRaw[1,:]<(IntFieldRaw[1,(Fixed_zLocations_IntField[i])]+(DipoleSeperation/2)))&((IntFieldRaw[2,(Fixed_zLocations_IntField[i])]-(DipoleSeperation/2))<IncBeamRaw[2,:])&(IncBeamRaw[2,:]<(IntFieldRaw[2,(Fixed_zLocations_IntField[i])]+(DipoleSeperation/2))))
    axis_array[i,0]=IncBeamRaw[0,IncBeamLineNumber]
    axis_array[i,1]=IncBeamRaw[1,IncBeamLineNumber]
    axis_array[i,2]=IncBeamRaw[2,IncBeamLineNumber]
    EField_array[i,0]=(IncBeamRaw[4,IncBeamLineNumber]+IntFieldRaw[4,(Fixed_zLocations_IntField[i])])+(1j*(IncBeamRaw[5,IncBeamLineNumber]+IntFieldRaw[5,(Fixed_zLocations_IntField[i])]))
    EField_array[i,1]=(IncBeamRaw[6,IncBeamLineNumber]+IntFieldRaw[6,(Fixed_zLocations_IntField[i])])+(1j*(IncBeamRaw[7,IncBeamLineNumber]+IntFieldRaw[7,(Fixed_zLocations_IntField[i])]))
    EField_array[i,2]=(IncBeamRaw[8,IncBeamLineNumber]+IntFieldRaw[8,(Fixed_zLocations_IntField[i])])+(1j*(IncBeamRaw[9,IncBeamLineNumber]+IntFieldRaw[9,(Fixed_zLocations_IntField[i])]))


plt.figure(1)

plt.title('Real Electric Field')

plt.quiver(axis_array[:,0], axis_array[:,1], np.real(EField_array[:,0]), np.real(EField_array[:,1]))

plt.figure(2)

plt.title('Imaginary Electric Field')

plt.quiver(axis_array[:,0], axis_array[:,1], np.imag(EField_array[:,0]), np.imag(EField_array[:,1]))

plt.show()