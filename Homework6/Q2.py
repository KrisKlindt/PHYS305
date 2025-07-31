# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 13:06:37 2024

@author: krist
"""
import numpy as np

def GaussJordanElimination(A,b): 
    # peforms Gauss-Jordan elimination to solve A x = b
    
    # Assuming A is numpy array
    # determine the shape of A
    L,W = A.shape

    # define matrix B so that A does not get overwritten
    B = np.array(A,dtype='float')
    
    # initialize array for solution x
    x = np.array(b,dtype='float')
    
    vals = np.array([k for k in range(W)])

    # begin iteration
    for j in range(W):
        # check if diagonal element is zero, pivot if it is
        if B[j,j] == 0:
            # in the jth column find and number points that are nonzero
            nonzero = [k for k,e in enumerate(B[j:,j]) if e != 0]
            # grab the number of the first nonzero element
            val = nonzero[0]
            # extract the current row and first row with nonzero diagonal
            b1 = B[j+val,:].copy()
            b2 = B[j,:].copy()
            # swap those rows
            B[j,:] = b1
            B[j+val,:] = b2
        
            # get the current row and first row with nonzero diagonal
            c1 = x[j+val].copy()
            c2 = x[j].copy()
            # swap those rows
            x[j] = c1
            x[j+val] = c2
            
        #divide the jth row of B by the diagonal element, likewise for q
        norm = B[j,j]
        B[j,:] = B[j,:]/norm
        x[j] = x[j]/norm

        # for all i not equal to j subtract correct multiple of jth row to
        # get zero in the jth column of the ith row
        for i in vals[(vals!=j)&(B[:,j]!=0)]:
            norm = B[i,j]
            B[i,:] = B[i,:] - norm*B[j,:]
            x[i] = x[i] - norm*x[j]
            
    return x

Equations = np.array([[-1,1,0,0,0],
                      [1,-3,1,0,0],
                      [0,1,-3,1,0],
                      [0,0,1,-3,1],
                      [0,0,0,-1,1]])

Solutions = np.array([-1,0,0,0,1])

Variables = GaussJordanElimination(Equations, Solutions)
print(Variables)