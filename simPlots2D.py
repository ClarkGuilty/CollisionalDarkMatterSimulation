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
rcParams.update({'font.size': 10})
plt.rcParams['image.cmap'] = 'plasma'
dpiT = 200
fsize = 16

def fastShow(image, title="none"):
    plt.figure()
    plt.imshow(image)
    plt.colorbar()
    plt.title(title)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

constantes = np.loadtxt("datFiles/constants.dat", usecols = 1)
Nt = int(constantes[-2])
TAU = int(constantes[-1])

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
fsize=16
interval = 1
dpII = 200
velUnit = 1183 #m/s
estUnit = 50 #kpc
potUnit = 1400318153625 #J/kg
acceUnit = 9.0761782e-13 #km/s²
dt = 0.5
plt.figure()
for i in range(0,Nt,interval):
    dens = np.loadtxt(giveDens(i,0)).T
    plt.imshow(dens,extent=[-1,1,-1,1],aspect = 'auto')
    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]",fontsize=fsize)
    plt.yticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
    plt.ylabel("Position [kpc]",fontsize=fsize)
#    plt.clim(0,1.2)
    cbar.set_label("Density [$M_{\odot}$ / kpc$^2$]",fontsize=fsize)
    #plt.title("Density $\\tau =$ {:d}".format(TAU),fontsize=fsize)
    plt.title("Density $t =$ {:.2f} ut".format(i*dt),fontsize=fsize)
    #plt.title('Densidad t = %d' %(i))
    plt.savefig(giveDens(i,1),dpi=dpi)
    plt.clf()

#
for i in range(0,Nt,interval):
    phasex = np.loadtxt(giveGridX(i,0)).T
    plt.imshow(phasex,extent=[-1,1,-1,1])
    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]",fontsize=fsize)
    plt.yticks(plt.xticks()[0], [str(t*velUnit) for t in plt.xticks()[0]])
    plt.ylabel("Velocity [km/s]",fontsize=fsize)
    plt.title("Phase Space $t =$ {:.2f} ut".format(i*dt),fontsize=fsize)
    cbar.set_label("Phase space density [$M_{\odot}$ / (kpc  $\\frac{km}{s}$)$^2$]",fontsize=fsize-2)
    plt.savefig(giveGridX(i,1),dpi=dpi)
    plt.clf()


#
#for i in range(0,Nt,interval):
#    phasey = np.loadtxt(giveGridY(i,0)).T
#    plt.imshow(phasey,extent=[-1,1,-1,1])
#    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
#    plt.xlabel("Position [kpc]",fontsize=fsize)
#    plt.yticks(plt.xticks()[0], [str(t*velUnit) for t in plt.xticks()[0]])
#    plt.ylabel("Velocity [km/s]",fontsize=fsize)
#    plt.title("Phase space $t =$ {:.2f} ut".format(i*dt),fontsize=fsize)
#    #plt.clim(0,1e-4)
#    cbar.set_label("Density [$M_{\odot}$ / kpc$^2$]",fontsize=fsize)
#    plt.savefig(giveGridY(i,1),dpi=dpi)
#    plt.clf()



























