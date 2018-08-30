#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 15:55:35 2018

@author: Javier Alejandro Acevedo Barroso
Script de Python para la visualización de la simulación en 2D.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

constantes = np.loadtxt("constants.dat", usecols = 1)
TAU = int(constantes[8])

#x = np.linspace(constantes[0], constantes[1], int(constantes[4]))  

densidad = np.loadtxt('./datFiles/density0.dat')
densidad1 = np.loadtxt('./datFiles/potential0.dat')

plt.imshow(densidad)

plt.figure()
plt.imshow(densidad1)

#def fmt(x, pos):
#    a, b = '{:.1e}'.format(x).split('e')
#    b = int(b)
#    return r'${} \times 10^{{{}}}$'.format(a, b)
#
#def giveDens(i,n):
#    cosa = ['a','b']
#    cosa[0] = './datFiles/density{:d}.dat'.format(i)
#    cosa[1] = './images/density{:d}.png'.format(i)
#    return cosa[n]
#
#def giveGridX(i,n):
#    cosa = ['a','b']
#    cosa[0] = './datFiles/gridx{:d}.dat'.format(i)
#    cosa[1] = './images/gridx{:d}.png'.format(i)
#    return cosa[n]
#
#def giveGridY(i,n):
#    cosa = ['a','b']
#    cosa[0] = './datFiles/gridy{:d}.dat'.format(i)
#    cosa[1] = './images/gridy{:d}.png'.format(i)
#    return cosa[n]
#
#
#dpi = 300
#
#for i in range(15):
#    dens = np.loadtxt(giveDens(i,0)).T
#    h0 = plt.figure()
#    plt.imshow(dens)
#    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#    plt.title('Densidad t = %d' %(i))
#    cbar.set_clim(0,4e-2)
#    plt.savefig(giveDens(i,1),dpi=dpi)
#
#
#for i in range(15):
#    phasex = np.loadtxt(giveGridX(i,0)).T
#    h0 = plt.figure()
#    plt.imshow(phasex)
#    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#    plt.title('phaseX t = %d' %(i))
#    plt.savefig(giveGridX(i,1),dpi=dpi)
#
#for i in range(15):
#    phasey = np.loadtxt(giveGridY(i,0)).T
#    h0 = plt.figure()
#    plt.imshow(phasey)
#    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#    plt.clim(0,1e-4)
#    plt.title('phaseY t = %d' %(i))
#    plt.savefig(giveGridY(i,1),dpi=dpi)
#
#
#
#
#



























