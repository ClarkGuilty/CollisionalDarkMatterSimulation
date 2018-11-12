#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 17:30:52 2018

@author: Javier Alejandro Acevedo Barroso
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams

fsize = 16
rcParams.update({'figure.autolayout': True})
#plt.rcParams['image.cmap'] = 'PuBu'
#plt.rcParams['image.cmap'] = 'YlGnBu'
rcParams.update({'font.size': 11})
plt.rcParams['image.cmap'] = 'plasma'

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

dat = np.loadtxt("./datFiles/grid0.dat").T
#density = np.loadtxt("density.dat")


grid0 = np.loadtxt("./datFiles/grid{:d}.dat".format(0)).T
dpII = 300

#        
#figu = plt.gcf()
figure = plt.figure(figsize=(7,5))
#figu.set_size_inches(18.5, 10.5)
#figu.set_dpi(300)


###########################################          1D ##############################


#constantes = np.loadtxt("constants.dat", usecols = 1)
#TAU = int(constantes[8])
#x = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
#velUnit = 1183 #m/s
#estUnit = 50 #kpc
#potUnit = 1400318153625 #J/kg
#acceUnit = 9.0761782e-13 #km/s²

#plt.imshow(dat, extent=[constantes[0],constantes[1],constantes[2],constantes[3]], interpolation='nearest', aspect='auto') #Es mucho más rápido imshow	
#plt.yticks(plt.yticks()[0], [str(np.round(t*velUnit)) for t in plt.yticks()[0]]) 
#plt.ylabel("Velocity [km/s]",fontsize=fsize)
#plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
#plt.xlabel("Position [kpc]", fontsize=fsize)
#plt.title("Phase Space Initialization",fontsize=fsize)        
#plt.clim(0,1e5)
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#cbar.set_label("Mass density [$M_{\odot}$ / kpc  $\\frac{km}{s}$]",fontsize=fsize+1)
#plt.savefig("1dInitPS.png", dpi = dpII)
#plt.clf()
#
#
#dens = np.loadtxt("./datFiles/density{:d}.dat".format(0))    
#plt.plot(x,dens)
#plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
#plt.xlabel("Position [kpc]",fontsize=fsize)
#plt.ylabel("Linear Density [$M_{\odot}$ / kpc]",fontsize=fsize)
#plt.title("Density Initialization",fontsize=fsize)        
#plt.ylim(-0.75e9,5.3e10)
#plt.xlim(-1.05,1.05)
#
#plt.savefig("1dInitDens.png", dpi = dpII)
#plt.clf()
#
# 
#potential = np.loadtxt("./datFiles/potential{:d}.dat".format(0))
#plt.plot(x,potential)
#plt.ylabel("Potential [J /kg]",fontsize=fsize)
#plt.title("Potential at t=0".format(TAU),fontsize=fsize)
##plt.ylim(-6.6e10,-5.8e10)
#plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
##plt.yticks(plt.yticks()[0], [fmt(np.round(t*potUnit),1) for t in plt.yticks()[0]]) 
#plt.xlabel("Position [kpc]",fontsize=fsize)
#plt.xlim(-1.05,1.05)
#plt.savefig("1dInitPot.png", dpi = dpII)
#plt.clf()
#
#
#acce = np.loadtxt("./datFiles/acce{:d}.dat".format(0))
#plt.plot(x,acce)
#plt.ylabel("Acceleration [kpc / Gy$^2$]",fontsize=fsize)
#plt.title("Acceleration at t=0".format(TAU),fontsize=fsize)
##plt.yticks(plt.yticks()[0], [str(t*2754463327) for t in plt.yticks()[0]])
#plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
##plt.yticks(plt.yticks()[0], [fmt(t*acceUnit,1) for t in plt.yticks()[0]]) 
#plt.ylim(np.min(acce)*1.1,np.max(acce)*1.1)
#
#plt.xlabel("Position [kpc]",fontsize=fsize)
#plt.xlim(-1.05,1.05)
#plt.savefig("1dInitAcce.png", dpi = dpII)
#plt.clf()

####################################################################


################################################################### 2D ########################


velUnit = 1183 #m/s
estUnit = 50 #kpc
potUnit = 1400318153625 #J/kg
acceUnit = 9.0761782e-13 #km/s²


dens = np.loadtxt("./miniCluster/2D/density0.dat").T
plt.imshow(dens,extent=[-1,1,-1,1])
cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
plt.xlabel("Position [kpc]",fontsize=fsize)
plt.yticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
plt.ylabel("Position [kpc]",fontsize=fsize)
cbar.set_label("Density [$M_{\odot}$ / kpc$^2$]",fontsize=fsize)
plt.title("2D Density Initialization",fontsize=fsize)
plt.savefig("2dInitDens.png", dpi = dpII)
plt.clf()

phasex = np.loadtxt('./miniCluster/2D/gridx0.dat').T
plt.imshow(phasex,extent=[-1,1,-1,1])
cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
plt.xlabel("Position [kpc]",fontsize=fsize)
plt.yticks(plt.xticks()[0], [str(t*velUnit) for t in plt.xticks()[0]])
plt.ylabel("Velocity [km/s]",fontsize=fsize)
plt.title("Phase Space initialization cut at y=0, Vy = 0",fontsize=fsize)
cbar.set_label("Phase Space Density [$M_{\odot}$ / (kpc  $\\frac{km}{s}$)$^2$]",fontsize=fsize)
plt.savefig('2dInitPhase.png',dpi=dpII)
plt.clf()














