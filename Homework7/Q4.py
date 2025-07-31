# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:39:36 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt

File = open('Problem4_data.txt','r')
Data = []
for line in File:
    Data.append(line.split())
    
#print(Data)

Data = np.array(Data, 'float')
N,d = Data.shape

fFFT = np.fft.fft(Data[:,1])

#plt.plot(fFFT)

freq = [2*np.pi*i/(Data[-1,0] - Data[0,0]) for i in range(N)]

#plt.plot(freq[0:round(N/2)], fFFT[0:round(N/2)])

Pow = np.conj(fFFT)*fFFT

#plt.plot(freq[0:round(N/2)], Pow[0:round(N/2)], 'b')

freq = np.array(freq)

Mask = (freq<10)

plt.plot(freq[Mask],Pow[Mask], 'b')