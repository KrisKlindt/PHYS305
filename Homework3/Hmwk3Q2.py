# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 17:30:15 2024

@author: krist
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np

def Gforce(m1, r1, m2, r2):
    G = 1
    r = np.sqrt(np.sum((r1 - r2) ** 2))
    F = -G * m1 * m2 / r**3 * (r2 - r1)
    return F

def Question2(TotalTime,dt,Skip):
# Computes the dynamics of the sun, Earth, Moon, and Jupiter
# Assumes orbits are circular

# define parameters
    mS = 1.0       # Mass of the Sun in solar masses
    mE = 3e-6    # Mass of Earth in solar masses
    mJ = 9.5e-4 # Mass of Jupiter in solar masses
    mM = 3.7e-8   # Mass of the Moon in solar masses
    
# define step size
    Steps = round(TotalTime/dt/Skip)

# initialize positions and velocities
    rS = np.zeros((Steps,2), dtype = np.float64)
    rE = np.zeros((Steps,2), dtype = np.float64)
    rM = np.zeros((Steps,2), dtype = np.float64)
    rJ = np.zeros((Steps,2), dtype = np.float64)
    
    vS = np.zeros((Steps,2), dtype = np.float64)
    vE = np.zeros((Steps,2), dtype = np.float64)
    vM = np.zeros((Steps,2), dtype = np.float64)
    vJ = np.zeros((Steps,2), dtype = np.float64)

    Time = np.zeros((Steps,1),dtype = np.float64)
    
    rS[0,0] = 0 
    rE[0,0] = 1 # in AU
    #for the Moon, will act as if Earth is (0,0) point
    rM[0,0] = 0.0026 +1# in AU
    rJ[0,0] = 5.2 # in AU
    
    vS[0,1] = 0
    vE[0,1] = 1 # arbitrary units
    vM[0,1] = 0.034 # scaled to Earth's velocity around sun
    vJ[0,1] = 0.436 # scaled to Earth's velocity around sun
    
    rSmid = np.zeros((1,2),dtype = np.float64)
    rEmid = np.zeros((1,2),dtype = np.float64)
    rMmid = np.zeros((1,2),dtype = np.float64)
    rJmid = np.zeros((1,2),dtype = np.float64)
    
    vSmid = np.zeros((1,2),dtype = np.float64)
    vEmid = np.zeros((1,2),dtype = np.float64)
    vMmid = np.zeros((1,2),dtype = np.float64)
    vJmid = np.zeros((1,2),dtype = np.float64)
    
    # movie stuff
    metadata = dict(title='Question2', artist='Matplotlib')
    writer = FFMpegWriter(fps=15,metadata= metadata)
    
    fig1 = plt.figure()
    l1, = plt.plot([], [],'bo')
    l2, = plt.plot([], [],'ro')
    l3, = plt.plot([], [],'go')
    l4, = plt.plot([], [],'yo')
    
    plt.xlim(-6, 6)
    plt.ylim(-6, 6)
    
    with writer.saving(fig1, "Question2.mp4", 100):
        
        for i in range(1,Steps):
        # start time stepping
            rS[i,:] = rS[i-1,:]
            vS[i,:] = vS[i-1,:]
            
            rE[i,:] = rE[i-1,:]
            vE[i,:] = vE[i-1,:]
            
            rM[i,:] = rM[i-1,:] + rE[i-1,:]
            vM[i,:] = vM[i-1,:]
            
            rJ[i,:] = rJ[i-1,:]
            vJ[i,:] = vJ[i-1,:]
            
            Time[i] = Time[i-1]
            
            for j in range(Skip):
                # compute r's at half time step
                rSmid[0,:] = rS[i,:] + dt*vS[i,:]/2
                rEmid[0,:] = rE[i,:] + dt*vE[i,:]/2
                rMmid[0,:] = rM[i,:] + dt*vM[i,:]/2 
                rJmid[0,:] = rJ[i,:] + dt*vJ[i,:]/2
                
                # define the force at the half time step\
                Fsj = Gforce(mS,rSmid[0,:],mJ,rJmid[0,:])# between jupiter and sun
                Fse = Gforce(mS,rSmid[0,:],mE,rEmid[0,:]) # between sun and earth
                Fem = Gforce(mE,rEmid[0,:],mM,rMmid[0,:]) # between earth and moon
                
                # time step the velocities a full time step by midpoint
                vSmid[0,:] = 0.5*vS[i,:]
                vEmid[0,:] = 0.5*vE[i,:]
                vMmid[0,:] = 0.5*vM[i,:]
                vJmid[0,:] = 0.5*vJ[i,:]
                
                vS[i,:] = vS[i,:] + dt*Fsj/mS
                vJ[i,:] = vJ[i,:] + dt*Fsj/mJ
                vE[i,:] = vE[i,:] + dt*Fse/mE
                vM[i,:] = vM[i,:] + dt*Fem/mM
                
                vSmid[0,:] += 0.5*vS[i,:]
                vEmid[0,:] += 0.5*vE[i,:]
                vMmid[0,:] += 0.5*vM[i,:]
                vJmid[0,:] += 0.5*vJ[i,:]
                
                # time step the positions a full time step by midpoint
                rS[i,:] = rS[i,:] + dt*vSmid[0,:]
                rE[i,:] = rE[i,:] + dt*vEmid[0,:]
                rM[i,:] = rM[i,:] + dt*vMmid[0,:]
                rJ[i,:] = rJ[i,:] + dt*vJmid[0,:]
                
            l1.set_data(rS[i,0],rS[i,1])
            l2.set_data(rE[i,0],rE[i,1])
            l3.set_data(rM[i,0],rM[i,1])
            l4.set_data(rJ[i,0],rJ[i,1])
            writer.grab_frame()
            plt.pause(0.02)
            
    return rS,rJ,Time

Question2(100, 0.1, 10)