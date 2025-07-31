# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:52:23 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt
    
def CalcRMag(r):
    # Takes the given R and calculates its magnitude
    return (np.sqrt(r[0]**2 + r[1]**2))

def Q1(m, ceta, alpha, beta, v0, dt, TotalTime):
    
    # define step size
    Steps = round(TotalTime/dt)
    
    v1 = np.zeros((Steps,2),'float')
    d1 = np.zeros((Steps,2),'float')
    r1 = np.zeros((Steps,2),'float')
    v2 = np.zeros((Steps,2),'float')
    d2 = np.zeros((Steps,2),'float')
    r2 = np.zeros((Steps,2),'float')

    Time = np.zeros((Steps,1),'float')
    
    # initial conditions
    v1[0,0] = v0
    d1[0,1] = 1 # moving in y direction
    r1[0,0] = 1
    v2[0,0] = v0
    d2[1,0] = 1 # moving in x direction
    r2[0,0] = 0
    
    
    for i in range(1,Steps):
        r = CalcRMag(abs(r1[i-1,:]-r2[i-1,:]))
        
        # Normalize the direction vectors to ensure they remain unit vectors
        d1[i,:] = d1[i,:] / np.linalg.norm(d1[i,:])
        d2[i,:] = d2[i,:] / np.linalg.norm(d2[i,:])
        
        # Use Forward Euler Method to evaluate v1, d1, v2, d2, r1, and r2 at next time step
        v1[i,:] = v1[i-1,:] + dt*((ceta*(v0*d1[i-1,:] - v1[i-1,:]) + (alpha/r)*(v2[i-1,:]-v1[i-1,:]))/m)
        a = np.cross(d1[i-1,:],d2[i-1,:])
        d1[i,:] = (beta/r**2)* a*d2[i-1,:]
        
        v2[i,:] = v2[i-1,:] + dt*((ceta*(v0*d2[i-1,:] - v2[i-1,:]) + (alpha/r)*(v1[i-1,:]-v2[i-1,:]))/m)
        b = np.cross(d2[i-1,:],d1[i-1,:])
        d2[i,:] = (beta/r**2)* b*d1[i-1,:]
        
        
        r1[i,:] = r1[i-1,:] + dt*v1[i-1,:]
        r2[i,:] = r2[i-1,:] + dt*v2[i-1,:]
        
        Time[i,0] = Time[i-1,0] + dt
        
    #plt.plot(r1[:,0], r1[:,1], label = "Particle 1")
    #plt.plot(r2[:,0], r2[:,1], label = "Particle 2")
    plt.plot(Time[:,0], v1[:,0], label = "v1x vs time")
    plt.plot(Time[:,0], v1[:,1], label = "v1y vs time")
    plt.plot(Time[:,0], v2[:,0], label = "v2x vs time")
    plt.plot(Time[:,0], v2[:,1], label = "v3x vs time")
    plt.legend()
        
        
Q1(1, 1, 1, 1, 0.5, 0.01, 100)
        