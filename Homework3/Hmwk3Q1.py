# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:21:42 2024

@author: krist
"""

# F = ma = -kR + q/r^3 R
# R is vector displacement with magnitude r, k/m = 1
# initial q=0, orbit circular

# m dv/dt = -kR
# dR/dt = v

# a = v**2/r ghvgh

import matplotlib.pyplot as plt
import numpy as np

def force(x, y, q):
    r = np.sqrt(x**2 + y**2)
    fx = -x + q * x / r**3
    fy = -y + q * y / r**3
    return fx, fy

def Spring(TotalTime, dt, Skip):
# Computes the dynamics of a spring

# define step size
    Steps = round(TotalTime/dt/Skip)

# initialize positions and velocities
# matrix with Steps rows and 2 columns
    r1 = np.zeros((Steps,2),'float')
    v1 = np.zeros((Steps,2),'float')

    Time = np.zeros((Steps,1),'float')
    
    # for this, q stands for q/m
    q = 2
    
    # initial conditions
    r1[0,0] = 1
    v1[0,1] = 1
    
    r1mid = np.zeros((1,2),'float')
    v1mid = np.zeros((1,2),'float')
    
    for i in range(1,Steps):
    # start time stepping
        r1[i,:] = r1[i-1,:]
        v1[i,:] = v1[i-1,:]
        Time[i] = Time[i-1]
            
        for j in range(Skip):
            # compute r1 and r2 at half time step
            r1mid[0,:] = r1[i,:] + dt*v1[i,:]/2
                
            # define the force at the half time step
            Fx, Fy = force(r1mid[0,0], r1mid[0,1],q)
                
            # time step the velocities a full time step by midpoint
            v1mid[0,:] = 0.5*v1[i,:]
            
            v1[i,0] = v1[i,0] + dt*Fx
            v1[i,1] = v1[i,1] + dt*Fy
                
            v1mid[0,:] += 0.5*v1[i,:]
                
            # time step the positions a full time step by midpoint
            r1[i,:] = r1[i,:] + dt*v1mid[0,:]
                    
    plt.plot(r1[:,0],r1[:,1],'b')
    return r1,v1,Time

Spring(30, 0.01, 10)
