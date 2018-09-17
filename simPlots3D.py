#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 12:32:13 2018

@author: Javier Alejandro Acevedo Barroso
Script de Python para la visualización de la simulación en 3D.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
dpiT = 500

def fastShow(image, title="none"):
    plt.clf()
    plt.imshow(image)
    cbar = plt.colorbar()
    plt.savefig(title+".png",dpi=dpiT)

    plt.title(title)

def fmt(x, pos):
    a, b = '{:.4e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

constantes = np.loadtxt("./datFiles/constants.dat", usecols = 1)
tau = int(constantes[-2])


densidadXY0 = np.loadtxt('./datFiles/densXY0.dat').T
densidadYZ0  = np.loadtxt('./datFiles/densYZ0.dat').T
densidadXZ0  = np.loadtxt('./datFiles/densXZ0.dat').T




#fastShow(densidadXY0, "XY")
#fastShow(densidadYZ0, "YZ")
#fastShow(densidadXZ0, "XZ")
plt.figure()
imagenes = 5
for i in range(0,tau,tau//imagenes):
    densidadXY = np.loadtxt("./datFiles/densXY{:d}.dat".format(i)).T
    potXY0 = np.loadtxt("./datFiles/pot0XY{:d}.dat".format(i)).T
    potXY1 = np.loadtxt("./datFiles/pot1XY{:d}.dat".format(i)).T
    diff = potXY0 - potXY1
    fastShow( diff, "XY{:d}".format(i) )
    