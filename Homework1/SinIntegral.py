# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 13:56:54 2024

@author: krist
"""

import math

def SinIntegral(LLimit,ULimit,Num):
# Computes the integral of sin(x) using the Left Hand, Right Hand, and
# Trapezoidal rule

# define step size
    Step = (ULimit-LLimit)/(Num-1)
    
# initialize Left Hand, Right Hand, and Trapezoidal rule Sum and Lower limit
    LHSum = 0
    RHSum = 0
    TrapSum = 0
    LowLimit = LLimit

# calculate Work done
    for i in range(Num-1):
        LHSum = LHSum + math.sin(LowLimit)*Step
        
        RHSum = RHSum + math.sin(LowLimit + Step)*Step
        
        TrapSum = TrapSum + (1/2)*(math.sin(LowLimit)+math.sin(LowLimit + Step))*Step
        
        LowLimit = LowLimit + Step
    
    return (LHSum, RHSum, TrapSum)

# set Num in increments of 10, up to 100
# each value in Num is used to determine step size of integral
Num = [10 + 10*i for i in range(10)]

# intialize lists to hold values for error of each method
LHErrorList = [0 for i in range(10)]
RHErrorList = [0 for i in range(10)]
TrapErrorList = [0 for i in range(10)]

# ask user to define lower and upper limits of integration
LLim = float(input("Enter lower limit of integration: "))
ULim = float(input("Enter upper limit of integration: "))

for i in range(10):
    
    LHValue, RHValue, TrapValue = SinIntegral(LLim, ULim, Num[i])
    
# calculate error for each step size, add to corresponding list
    LHErrorList[i] = abs(LHValue-(-math.cos(ULim)+math.cos(LLim)))
    RHErrorList[i] = abs(RHValue-(-math.cos(ULim)+math.cos(LLim)))
    TrapErrorList[i] = abs(TrapValue-(-math.cos(ULim)+math.cos(LLim)))
    
lnN = [ math.log(Num[i]) for i in range(10) ]
lnLHE = [ math.log(LHErrorList[i]) for i in range(10) ]
lnRHE = [ math.log(RHErrorList[i]) for i in range(10) ]
lnTE = [ math.log(TrapErrorList[i]) for i in range(10) ]