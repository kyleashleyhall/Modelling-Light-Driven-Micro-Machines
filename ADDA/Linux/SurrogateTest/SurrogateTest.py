# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 12:03:16 2018

@author: Marius & Kyle
"""

import numpy as np
from da.p7core import gtapprox

def DipSep(Singleaxis):
    dx = np.zeros([len(Singleaxis)-1])
    for i in range(len(Singleaxis)-1):
        dx[i] = Singleaxis[i+1] - Singleaxis[i]
        if (dx[i]==0):
            dx[i]=np.inf
    return min(np.absolute(dx))

DipPol=np.transpose(np.loadtxt('DipPol-Y', skiprows=1))
IntField=np.transpose(np.loadtxt('IntField-Y', skiprows=1))
IncBeam=np.transpose(np.loadtxt('IncBeam-Y', skiprows=1))
DipoleSep=DipSep(DipPol[0,:])

#To be integrated into Forces function:

TotalField=np.zeros([9,len(DipPol[0])])


for i in range(len(DipPol[0])):

    linePositionCurrent_IncBeam=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IncBeam[0,:])&(IncBeam[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IncBeam[1,:])&(IncBeam[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IncBeam[2,:])&(IncBeam[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
    linePositionCurrent_IntField=np.argwhere(((DipPol[0,i]-((1*DipoleSep)/2))<IntField[0,:])&(IntField[0,:]<(DipPol[0,i]+((1*DipoleSep)/2)))&((DipPol[1,i]-((1*DipoleSep)/2))<IntField[1,:])&(IntField[1,:]<(DipPol[1,i]+((1*DipoleSep)/2)))&((DipPol[2,i]-((1*DipoleSep)/2))<IntField[2,:])&(IntField[2,:]<(DipPol[2,i]+((1*DipoleSep)/2))))
    TotalField[0,i]=DipPol[0,i]
    TotalField[1,i]=DipPol[1,i]
    TotalField[2,i]=DipPol[2,i]
    TotalField[3,i]=IncBeam[4,linePositionCurrent_IncBeam]+IntField[4,linePositionCurrent_IntField]
    TotalField[4,i]=IncBeam[5,linePositionCurrent_IncBeam]+IntField[5,linePositionCurrent_IntField]
    TotalField[5,i]=IncBeam[6,linePositionCurrent_IncBeam]+IntField[6,linePositionCurrent_IntField]
    TotalField[6,i]=IncBeam[7,linePositionCurrent_IncBeam]+IntField[7,linePositionCurrent_IntField]
    TotalField[7,i]=IncBeam[8,linePositionCurrent_IncBeam]+IntField[8,linePositionCurrent_IntField]
    TotalField[8,i]=IncBeam[9,linePositionCurrent_IncBeam]+IntField[9,linePositionCurrent_IntField]

xSample=np.transpose(TotalField[0:3,:])
ySample=np.transpose(TotalField[3,:])

print(xSample)
print(ySample)

builder = gtapprox.Builder()
options = {
'GTApprox/AccuracyEvaluation': 'on',
'GTApprox/Technique': 'GP',
'GTApprox/LogLevel': 'Info'
}
builder.options.set(options)
model = builder.build(xSample, ySample)

model.save('approxModel.gta')