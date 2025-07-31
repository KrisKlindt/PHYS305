# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 18:03:24 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt

def Brownian():
    # Physical constants
    kBT =  4*10**(-14)  # Thermal Energy
    eta = 0.01  # Water viscosity
    d = 2e-6  # Diameter of the particle in meters
    r = d / 2  # Radius of the particle in meters
    m = (4/3) * np.pi * r**3  # Mass of the particle
    zeta = 6 * np.pi * eta * r  # Drag coefficient
    dt = 1e-9  # Time step in s
    c = zeta*np.sqrt(kBT/m) # magnitude of stochastic force
    
    TotalTime = 1e-5  # Total simulation time (s)
    Steps = int(TotalTime / dt)  # Number of time steps
    
    # initialize lists for coordinates and velocities
    r = np.zeros((2,Steps))
    v = np.zeros((2,Steps))
    Time = [ dt*i for i in range(Steps) ]
    MSD = np.zeros((1,Steps))
    MSDlog = np.zeros((1,Steps))
    Timelog = [np.log(dt*i) for i in range(Steps)]

    # Simulate Brownian motion using the midpoint method
    for i in range(1, Steps):
        F = c* np.random.randn(2)
        
        v[:, i] = ((1 - zeta * dt /2/m) * v[:, i-1] + F * dt / m)/(1 + zeta * dt/2/m)
        
        vMid = 0.5 * (v[:,i] + v[:,i-1])

        r[:, i] = r[:, i-1] + vMid * dt
        
        MSD[0,i] = (MSD[0,i-1] + r[0, i]**2 + r[1, i]**2)
        MSDlog[0,i] = np.log(MSD[0,i])

    plt.plot(Timelog, MSDlog[0,:], label = "Brownian Particle Mean Square Displacement")
    plt.legend()
    
Brownian()
        