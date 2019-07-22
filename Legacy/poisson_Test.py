#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 18:44:29 2019

@author: Javier Alejandro Acevedo Barroso
"""


import numpy as np
import matplotlib.pyplot as plt


densTheo = np.loadtxt("testDensity.dat")
potTheo = np.loadtxt("testPotential.dat")
potSolved = np.loadtxt("solvedPotential.dat")


#potTheo -= np.max(potTheo)
#potSolved -= np.min(potSolved)


dif = potTheo-potSolved
sumD = dif[potTheo!=0]/potTheo[potTheo!=0]/np.max(potTheo)
sumDT = sumD.sum()



plt.plot(potTheo)
plt.plot(potSolved)

plt.figure()
plt.plot(dif)

pasos = np.diff(potSolved*2/16584)  
pasos2 = np.diff(pasos*2/16584)

h=2/16584
f = potTheo
pasos2 = (f[2:-1] -2*f[1:-2] + f[0:-3])/(h**2)
plt.figure()
plt.plot(densTheo)
plt.plot(pasos2)