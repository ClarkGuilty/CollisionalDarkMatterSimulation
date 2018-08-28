#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 15:55:35 2018

@author: Javier Alejandro Acevedo Barroso
Script de Python para la visualización de la simulación en 2D.
"""
import numpy as np
import matplotlib.pyplot as plt

potTest0 = np.loadtxt('Density0.dat').T
potTest1 = np.loadtxt('Density1.dat').T

h0 = plt.figure()
plt.imshow(potTest0)
cbar = plt.colorbar()

h1 = plt.figure()
plt.imshow(potTest1)
cbar = plt.colorbar()
#
#
#diff = potTest1 - potTest0
#h2 = plt.figure()
#plt.imshow(diff)
#cbar = plt.colorbar()
#
#print(np.sum(np.sum(potTest0-potTest1)))

pxTest0 = np.loadtxt('px0.dat').T
pxTest1 = np.loadtxt('pxAft.dat').T


h0 = plt.figure()
plt.imshow(pxTest0)
cbar = plt.colorbar()

h1 = plt.figure()
plt.imshow(pxTest1)
cbar = plt.colorbar()