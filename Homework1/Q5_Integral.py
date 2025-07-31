# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 09:28:04 2024

@author: krist
"""

import math

# define step size
l = 20
dV = 0.1

r = math.floor(l/dV)

TrapSum = 0
# Cannot start at 0, since it would divide by zero, so start very close to zero instead
V = 0.000001

for i in range(r):
    
    TrapSum += (1/2)*(((1-math.exp(-V))/(V*(1 + V**2))) + ((1-math.exp(-(V+dV)))/((V+dV)*(1+(V+dV)**2))))*dV
    
    V = V + dV
    
print(TrapSum)