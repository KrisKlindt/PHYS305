# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 13:19:25 2024

@author: krist
"""

import math
from IdealGas import IdealGasWork

# define parameters using cgs units
N = 2.2*10**23
kBT = 4*10**(-14)
Vi = 1
Vf = 10
    
# set Num in increments of 10, up to 100
# each value in Num is used to determine step size of integral
Num = [10 + 10*i for i in range(10)]

# intialize lists to hold values for error of each method
LHErrorList = [0 for i in range(10)]
RHErrorList = [0 for i in range(10)]
TrapErrorList = [0 for i in range(10)]

for i in range(10):
    
    LHWork, RHWork, TrapWork = IdealGasWork(Vi,Vf,Num[i])
    
# calculate error for each step size, add to corresponding list
    LHErrorList[i] = abs(LHWork-(N*kBT)*math.log(Vf/Vi))
    RHErrorList[i] = abs(RHWork-(N*kBT)*math.log(Vf/Vi))
    TrapErrorList[i] = abs(TrapWork-(N*kBT)*math.log(Vf/Vi))
    
lnN = [ math.log(Num[i]) for i in range(10) ]
lnLHE = [ math.log(LHErrorList[i]) for i in range(10) ]
lnRHE = [ math.log(RHErrorList[i]) for i in range(10) ]
lnTE = [ math.log(TrapErrorList[i]) for i in range(10) ]
        