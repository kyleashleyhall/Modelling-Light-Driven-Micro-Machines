# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 17:29:24 2018

@author: Marius
"""

import numpy as np

def DipSep(Singleaxis):
    dx = np.zeros([len(Singleaxis)-1])
    for i in range(len(Singleaxis)-1):
        dx[i] = Singleaxis[i+1] - Singleaxis[i]
        if (dx[i]==0):
            dx[i]=np.inf
    return min(np.absolute(dx))
    
IntField = np.loadtxt('IntField-Y',skiprows=1)
IncBeam = np.loadtxt('IncBeam-Y',skiprows=1)
DipPol = np.loadtxt('DipPol-Y', skiprows=1)

#print(np.shape(IntField))


def Forces(Intfield, Incbeam, Dippol, Dipsep):
    Force = np.zeros((len(Intfield),7))
    for i in range(len(Force)):
        Indexvals = np.arange(0)
        Force[i,0:3] = Intfield[i,0:3]
        #Force in x
        for j in range(len(Force)):
            if Intfield[i,1] == Intfield[j,1] and Intfield[i,2] == Intfield[j,2] and (Intfield[j,0] - Dipsep) <= Intfield[i,0] <= (Intfield[j,0] + Dipsep):
                    Indexvals = np.append(Indexvals, j)                
        for p in range(len(Indexvals)):
            if Intfield[i,0] == Intfield[Indexvals[p],0]:
                if len(Indexvals) == 2:
                    dx = (Intfield[Indexvals[-1]]+Incbeam[Indexvals[-1]]) - (Intfield[Indexvals[0]]+Incbeam[Indexvals[0]])
                    term1 = 2*(dx[4]+1j*dx[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dx[0]
                    term2 = 2*(dx[6]+1j*dx[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dx[0]
                    term3 = 2*(dx[8]+1j*dx[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dx[0]
                    Force[i,4] = 0.5*np.real(term1+term2+term3)
                
                elif len(Indexvals) == 3:
                    dx = (Intfield[Indexvals[-1]]+Incbeam[Indexvals[-1]]) - (Intfield[Indexvals[0]]+Incbeam[Indexvals[0]])
                    term1 = (dx[4]+1j*dx[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dx[0]
                    term2 = (dx[6]+1j*dx[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dx[0]
                    term3 = (dx[8]+1j*dx[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dx[0]
                    Force[i,4] = 0.5*np.real(term1+term2+term3)
                    
                elif len(Indexvals) == 1:
                    dx = Intfield[Indexvals[0]]+Incbeam[Indexvals[0]]
                    term1 = 2*(dx[4]+1j*dx[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dx[0]
                    term2 = 2*(dx[6]+1j*dx[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dx[0]
                    term3 = 2*(dx[8]+1j*dx[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dx[0]
                    Force[i,4] = 0.5*np.real(term1+term2+term3)
        #Force in y         
        Indexvals = np.arange(0)
        for j in range(len(Force)):
            if Intfield[i,0] == Intfield[j,0] and Intfield[i,2] == Intfield[j,2] and (Intfield[j,1] - Dipsep) <= Intfield[i,1] <= (Intfield[j,1] + Dipsep):
                    Indexvals = np.append(Indexvals, j)                
        for p in range(len(Indexvals)):
            if Intfield[i,1] == Intfield[Indexvals[p],1]:
                if len(Indexvals) == 2:
                    dy = (Intfield[Indexvals[-1]]+Incbeam[Indexvals[-1]]) - (Intfield[Indexvals[0]]+Incbeam[Indexvals[0]])
                    term1 = 2*(dy[4]+1j*dy[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dy[1]
                    term2 = 2*(dy[6]+1j*dy[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dy[1]
                    term3 = 2*(dy[8]+1j*dy[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dy[1]
                    Force[i,5] = 0.5*np.real(term1+term2+term3)
                
                elif len(Indexvals) == 3:
                    dy = (Intfield[Indexvals[-1]]+Incbeam[Indexvals[-1]]) - (Intfield[Indexvals[0]]+Incbeam[Indexvals[0]])
                    term1 = (dy[4]+1j*dy[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dy[1]
                    term2 = (dy[6]+1j*dy[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dy[1]
                    term3 = (dy[8]+1j*dy[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dy[1]
                    Force[i,5] = 0.5*np.real(term1+term2+term3)
                    
                elif len(Indexvals) == 1:
                    dy = Intfield[Indexvals[0]]+Incbeam[Indexvals[0]]
                    term1 = 2*(dy[4]+1j*dy[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dy[1]
                    term2 = 2*(dy[6]+1j*dy[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dy[1]
                    term3 = 2*(dy[8]+1j*dy[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dy[1]
                    Force[i,5] = 0.5*np.real(term1+term2+term3)
        #Force in z
        Indexvals = np.arange(0)
        for j in range(len(Force)):
            if Intfield[i,1] == Intfield[j,1] and Intfield[i,0] == Intfield[j,0] and (Intfield[j,2] - Dipsep) <= Intfield[i,2] <= (Intfield[j,2] + Dipsep):
                    Indexvals = np.append(Indexvals, j)                
        for p in range(len(Indexvals)):
            if Intfield[i,2] == Intfield[Indexvals[p],2]: #THis condition finds the actual dipole labelled i in Indexvals
                if len(Indexvals) == 2:
                    dz = (Intfield[Indexvals[-1]]+Incbeam[Indexvals[-1]]) - (Intfield[Indexvals[0]]+Incbeam[Indexvals[0]])
                    term1 = 2*(dz[4]+1j*dz[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dz[2]
                    term2 = 2*(dz[6]+1j*dz[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dz[2]
                    term3 = 2*(dz[8]+1j*dz[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dz[2]
                    Force[i,6] = 0.5*np.real(term1+term2+term3)
                
                elif len(Indexvals) == 3:
                    dz = (Intfield[Indexvals[-1]]+Incbeam[Indexvals[-1]]) - (Intfield[Indexvals[0]]+Incbeam[Indexvals[0]])
                    term1 = (dz[4]+1j*dz[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dz[2]
                    term2 = (dz[6]+1j*dz[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dz[2]
                    term3 = (dz[8]+1j*dz[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dz[2]
                    Force[i,6] = 0.5*np.real(term1+term2+term3)
                    
                elif len(Indexvals) == 1:
                    dz = Intfield[Indexvals[0]]+Incbeam[Indexvals[0]]
                    term1 = 2*(dz[4]+1j*dz[5])*(Dippol[Indexvals[p],4]+1j*Dippol[Indexvals[p],5])/dz[2]
                    term2 = 2*(dz[6]+1j*dz[7])*(Dippol[Indexvals[p],6]+1j*Dippol[Indexvals[p],7])/dz[2]
                    term3 = 2*(dz[8]+1j*dz[9])*(Dippol[Indexvals[p],8]+1j*Dippol[Indexvals[p],9])/dz[2]
                    Force[i,6] = 0.5*np.real(term1+term2+term3)
        Force[i,3] = (Force[i,4]**2+Force[i,5]**2+Force[i,6]**2)
        print(i)
    return Force

C = Forces(IntField, IncBeam, DipPol, DipSep(IntField[:,0]))        
print(C, sum(C))
print(sum(C)[3]-sum(C)[-1]**2)
Rad = np.loadtxt('RadForce-Y', skiprows=1)
print(sum(Rad))
print(sum(Rad)[3]-sum(Rad[-1]**2))
    
    
    