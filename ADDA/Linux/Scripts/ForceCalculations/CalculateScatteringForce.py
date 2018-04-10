# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 14:27:42 2018

@author: Marius
"""

import scipy.constants as cons
import numpy as np


def svalue(wl, bw):
    return (wl/bw)/(2*np.pi)
    
def Intensity(power, s, permittivity, bw):
    return (16*power)/(np.sqrt(permittivity)*cons.c*(bw**2)*(1+s**2+(1.5*s**4)))

def ScatteringForce(m_medium, m_relative, wavelength, r, Intensity):
    Term1 = (8*cons.pi*m_medium*(((2*cons.pi)/wavelength)**4)*(r**6))/(3*cons.c)
    Term2 = Intensity*((((m_relative**2)-1)/((m_relative**2)+2))**2)
    return Term1*Term2
    
m_medium = 1.326
m_relative = 1.18595
Temperature = 20
wave = 1.064e-6
r = 1e-6
power = 5e-3
permittivity = (87.740-(0.40008*Temperature)+(9.398e-4*(Temperature**2))-(1.410e-6*(Temperature**3)))*(8.85e-12)

beamwidth=1e-6
s = svalue(wave, beamwidth)
print(s)
Intensity = Intensity(power, s, permittivity, beamwidth)
#Intensity = power/(np.pi*((beamwidth/2)**2))
print(Intensity)
print(ScatteringForce(m_medium, m_relative, wave, r, Intensity))