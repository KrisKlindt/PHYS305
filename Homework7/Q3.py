# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:09:18 2024

@author: krist
"""
import numpy as np
import matplotlib.pyplot as plt


def FitData(datafile):
    File = open(datafile,'r')
    Dat = []
    for line in File:
        Dat.append(line.split())
        
    #print(Data)
    Data = np.array(Dat, 'float')
    N,d = Data.shape
    
    lam = .01
    x = Data[:,0]
    dx = x[1] - x[0]

    M = np.zeros((N,N))
    b = Data[:,1].copy()
    b[0] = 0

    M[0,0] = -lam/dx
    M[0,1] = lam/dx

    for i in range(1,N-1):
        M[i,i-1] = lam/(dx**2)
        M[i,i] = -1 - 2*lam/(dx**2)
        M[i,i+1] = lam/(dx**2)
        
    M[-1,-2] = -lam/dx
    M[-1,-1] = lam/dx
    
    a = np.linalg.solve(M,b) * -1
    
    plt.plot(Data[:,0], Data[:,1], 'b')
    plt.plot(x, a, 'r')
    
    return a
    
a = FitData('Problem2_data.txt')