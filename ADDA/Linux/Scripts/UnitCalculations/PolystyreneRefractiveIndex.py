# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:26:28 2018

@author: Kyle
"""

def PolystyreneRefractiveIndex(Lambda,MediumsRefractiveIndex):
    
    RefracIndex=((((1.4435*((Lambda)**2)))/(((Lambda)**2)-0.020216))+1)**0.5
    RefracIndex/=MediumsRefractiveIndex
    
    return(RefracIndex)
    
print(PolystyreneRefractiveIndex(1.064,1.328))