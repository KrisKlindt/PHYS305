# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:07:44 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt

def FitDatatoCubic(datafile):
    File = open(datafile,'r')
    Dat = []
    for line in File:
        Dat.append(line.split())
        
    #print(Data)
    Data = np.array(Dat, 'float')
    N,d = Data.shape
    
    x = Data[:,0]
    
    M = np.zeros((4,4))
    M[:,0] = [N, np.sum(x), np.sum(x**2), np.sum(x**3)]
    M[:,1] = [np.sum(x), np.sum(x**2), np.sum(x**3), np.sum(x**4)]
    M[:,2] = [np.sum(x**2), np.sum(x**3), np.sum(x**4), np.sum(x**5)]
    M[:,3] = [np.sum(x**3), np.sum(x**4), np.sum(x**5), np.sum(x**6)]
    
    b = np.zeros((4,1))
    b[:,0] = [np.sum(Data[:,1]), np.sum(x*Data[:,1]), np.sum(x**2*Data[:,1]), np.sum(x**3*Data[:,1])]
    
    a = np.linalg.solve(M,b)
    
    y = a[0] + a[1]*x + a[2]*x**2 + a[3]*x**3
    
    plt.plot(Data[:,0], Data[:,1], 'b')
    plt.plot(x, y, 'r')
    
    return a
    
a = FitDatatoCubic('Problem2_data.txt')