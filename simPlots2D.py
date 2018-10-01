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


def fastShow(image, title="none"):
    plt.figure()
    plt.imshow(image)
    plt.colorbar()
    plt.title(title)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

constantes = np.loadtxt("constants.dat", usecols = 1)
Nt = int(constantes[-1])

#x = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
dpi = 200

#densidadTheo = np.loadtxt('./datFiles/density0.dat').T
#potTheo = np.loadtxt('./datFiles/potential0.dat').T
#potReal = np.loadtxt('./datFiles/potential1.dat').T


#plt.imshow(densidadTheo)
#plt.title("Densidad cumple la ecuación de poisson del potTeorico")
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#plt.savefig("densidadTeorica.png",dpi=dpi)

#plt.figure()
#plt.imshow(potTheo)
#plt.title("potTeorico V = cos(0.5*pi*x)cos(0.5*pi*y)")
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#plt.savefig("PotencialTeorico.png",dpi=dpi)

#plt.figure()
#plt.imshow(potReal)
#plt.title("potCalculado a partir de la densidad")
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#plt.savefig("potReal.png",dpi=dpi)


def darX(inx):
    return -1.0+2.0/128*inx

def darY(iny):
    return -1.0+2.0/128*iny

dista = (int) (0.2*Nt)
mina = (int) (Nt//2 ) - dista 
maxa = (int) ( Nt//2+dista) 


#diff = potTheo-potReal
#
#plt.figure()
#plt.imshow(diff)
#plt.title("potTeorico/potCalculado")
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#plt.savefig("potCalpotTeo.png", dpi=dpi)
#
#
#realx = np.loadtxt('./datFiles/realx.dat').T
#realy = np.loadtxt('./datFiles/realy.dat').T
#calcx = np.loadtxt('./datFiles/calcx.dat').T
#calcy = np.loadtxt('./datFiles/calcy.dat').T
#
#fastShow(realx, title="realx")
#fastShow(realy, title="realy")
#fastShow(calcx, title="accex Calc - Teorica")
#fastShow(calcy, title="accey Calc - Teorica")


#
    
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
def giveDens(i,n):
    cosa = ['a','b']
    cosa[0] = './datFiles/density{:d}.dat'.format(i)
    cosa[1] = './images/density{:d}.png'.format(i)
    return cosa[n]
#
def giveGridX(i,n):
    cosa = ['a','b']
    cosa[0] = './datFiles/gridx{:d}.dat'.format(i)
    cosa[1] = './images/gridx{:d}.png'.format(i)
    return cosa[n]
#
def giveGridY(i,n):
    cosa = ['a','b']
    cosa[0] = './datFiles/gridy{:d}.dat'.format(i)
    cosa[1] = './images/gridy{:d}.png'.format(i)
    return cosa[n]
#
#
#dpi = 300
#
interval = 1
for i in range(0,Nt,interval):
    dens = np.loadtxt(giveDens(i,0)).T
    h0 = plt.figure()
    plt.imshow(dens)
    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
    plt.title('Densidad t = %d' %(i))
    #cbar.set_clim(0,16)
    plt.savefig(giveDens(i,1),dpi=dpi)
#
#
for i in range(0,Nt,interval):
    phasex = np.loadtxt(giveGridX(i,0)).T
    h0 = plt.figure()
    plt.imshow(phasex)
    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
    plt.title('phaseX t = %d' %(i))
    plt.savefig(giveGridX(i,1),dpi=dpi)
#
for i in range(0,Nt,interval):
    phasey = np.loadtxt(giveGridY(i,0)).T
    h0 = plt.figure()
    plt.imshow(phasey)
    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
    #plt.clim(0,1e-4)
    plt.title('phaseY t = %d' %(i))
    plt.savefig(giveGridY(i,1),dpi=dpi)
#
#
#
#
#
#for i in range(128):
#    print(-1.0+i*2/127)



























