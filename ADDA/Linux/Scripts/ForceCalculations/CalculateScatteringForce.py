# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 14:27:42 2018

@author: Marius
"""

import scipy.constants as cons

def ScatteringForce(m_medium, m_relative, wavelength, r, Intensity):
    Term1 = (8*cons.pi*m_medium*(((2*cons.pi)/wavelength)**4)*(r**6))/(3*cons.c)
    Term2 = Intensity*((((m_relative**2)-1)/((m_relative**2)+2))**2)
    return Term1*Term2
    
m_medium = 1.326
m_relative = 1.18595
wave = 1.064e-6
r = 1e-6
Intensity = 1e7
print(ScatteringForce(m_medium, m_relative, wave, r, Intensity))