# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 21:20:23 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt

def Calch(dtheta, ai2, hi):
    k2 = dtheta * (ai2)
    
    return hi + k2
    
def CalcParams(ai, hi, theta, dtheta, z_i, m_i):
    distance_i = np.zeros_like(z_i)
    z0_minus_zi = z_i
    #x = h[i-1,0] * np.cos(theta)
    #z = h[i-1,0] * np.sin(theta)
    for c in range(len(z_i)):
        distance_i[c] = np.sqrt(hi ** 2 +
                               2.0 * z0_minus_zi[c] * hi * np.cos(theta)
                               + z0_minus_zi[c] ** 2)
        #print("distance " + str(distance_i[c]))
       
    psi = 1.0
    dpsi_dr = 0.0
    dpsi_dtheta = 0.0
    for b in range(len(m_i)):
        psi += 0.5 * m_i[b] / distance_i[b]
        dpsi_dr -= 0.5 * m_i[b] * (hi + z0_minus_zi[b] * np.cos(theta)) /\
            distance_i[b] ** 3
        dpsi_dtheta += 0.5 * m_i[b] * hi * z0_minus_zi[b] * np.sin(theta) /\
           distance_i[b] ** 3
       
    C2 = 1.0 / (1.0 + (ai / hi) ** 2)
       
    if (abs(theta) <= np.pi/100) or (abs(theta - np.pi) <= np.pi/100) \
        or (abs(2*np.pi - theta) <= np.pi/100):
        cot_theta_dhdtheta_C2 = 0.0
    else:
        cot_theta_dhdtheta_C2 = ai / (np.tan(theta) * C2)
           
    return psi, dpsi_dr, dpsi_dtheta, C2, cot_theta_dhdtheta_C2

def CalcA(ai, hi, theta, dtheta, z_i, m_i):
    psi, dpsi_dr, dpsi_dtheta, C2, cot_theta_dhdtheta_C2 = \
        CalcParams(ai, hi, theta, dtheta, z_i, m_i)
    
    k1 = dtheta * (2*hi - cot_theta_dhdtheta_C2 + \
                (4.0 * hi ** 2 / (psi * C2)) * \
                (dpsi_dr - dpsi_dtheta * ai / hi ** 2) + \
                3.0 * ai ** 2 / hi)
    
    psi, dpsi_dr, dpsi_dtheta, C2, cot_theta_dhdtheta_C2 = \
        CalcParams((ai + k1/2), hi, (theta + (dtheta/2)), dtheta, z_i, m_i)
        
    k2 = dtheta * (2*hi - cot_theta_dhdtheta_C2 + \
                (4.0 * hi ** 2 / (psi * C2)) * \
                (dpsi_dr - dpsi_dtheta * (ai + (k1/2)) / hi ** 2) + \
                3.0 * (ai + (k1/2)) ** 2 / hi)
    return ai + k2

def CalcHorizon(self, dtheta, FullCircle):
    """
    Parameters
    ----------
    dtheta : float
        This is the size of the change in theta for each step
    FullCircle : float
        This is what theta will go to to complete a full circle = 2pi

    Returns
    -------
    None.

    """
    # define step size
    Steps = round(FullCircle/dtheta)
    
    a = np.zeros((Steps+1,1),'float')
    h = np.zeros((Steps+1,2),'float')
    
    # initial conditions
    a[0,0] = 0
    h[0,0] = self.r0[0]
    theta = 0 + dtheta
    
    #print(self.r0[0])
    
    z_i = self.spacetime.z_positions
    m_i = self.spacetime.masses
    
    for i in range(1,Steps+1):
        a2 = CalcA(a[i-1,0], h[i-1,0], theta, dtheta, z_i, m_i)
        
        h2 = Calch(dtheta, a2, h[i-1,0])
                
        if (i == 10000 or i == 20000):
            a[i,0] = 0
            h[i,0] = self.r0[0]
        else:
            a[i,0] = a2
                
            #print(str(a[i,0]) + " a")
            h[i,0] = h2
            #print(str(h[i,0]) + " h " + str(i))
        h[i,1] = theta
            
        theta += dtheta
    
    self.H = h
    return None
        
        
CalcHorizon(np.pi/100, 2*np.pi)