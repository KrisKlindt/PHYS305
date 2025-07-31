# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 13:02:43 2024

@author: krist
"""

import math

def IdealGasWork(Vi,Vf,Num):
# Computes the work done by an ideal gas when it# expands 
#from Vi to Vf using the lefthand rule

# define parameters using cgs units
    N = 2.2*10**23
    kBT = 4*10**(-14)
    
# define step size
    dV = (Vf-Vi)/(Num-1)
    
# initialize Left Hand, Right Hand, and Trapezoidal rule work and Volume
    LHWork = 0
    RHWork = 0
    TrapWork = 0
    V = Vi

# calculate Work done
    for i in range(Num-1):
        LHWork = LHWork + (N*kBT/V)*dV
        
        RHWork = RHWork + (N*kBT/(V+dV))*dV
        
        TrapWork = TrapWork + (1/2)*((N*kBT/V)+(N*kBT/(V+dV)))*dV
        
        V = V + dV
    
    return (LHWork, RHWork, TrapWork)