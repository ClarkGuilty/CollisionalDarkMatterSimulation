# -*- coding: utf-8 -*-
"""


@author: Javier Alejandro Acevedo Barroso
Script de Python para la visualización de la simulación.

"""
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib.ticker as ticker
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

#plt.rcParams['image.cmap'] = 'PuBu'
#plt.rcParams['image.cmap'] = 'YlGnBu'
plt.rcParams['image.cmap'] = 'plasma'
rcParams.update({'font.size': 11})
plt.rcParams['image.cmap'] = 'plasma'
fsize = 16

BULLET = -147
JEANS = -137
GAUSS = -127


dat = np.loadtxt("./datFiles/grid0.dat").T
#density = np.loadtxt("density.dat")
constantes = np.loadtxt("./datFiles/constants.dat", usecols = 1)
dt = constantes[9]
TAU = int(constantes[8])
#inF = np.loadtxt("inF.dat")
#outF = np.loadtxt("outF0.dat")
#outF1 = np.loadtxt("outF1.dat")
#oI = np.loadtxt("oI.dat")
#oR = np.loadtxt("oR.dat")
#acce = np.loadtxt("acce.dat")

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

x = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
#        
figu = plt.gcf()
#figu.set_size_inches(18.5, 10.5)
#figu.set_dpi(300)
dpII = 200
velUnit = 621 #m/s
estUnit = 35 #kpc
potUnit = 385962691092 #J/kg
acceUnit = 3.5737451e-13 #km/s²


for i in range(int(constantes[6])):
    dat = np.loadtxt("./datFiles/grid{:d}.dat".format(i)).T
    #dat = dat#/np.max(dat)/7
    plt.imshow(dat, extent=[constantes[0],constantes[1],constantes[2],constantes[3]], aspect='auto') #Es mucho más rápido imshow
#    plt.contourf(np.flip(dat,axis=0), extent=[constantes[0],constantes[1],constantes[2],constantes[3]], levels = 8) #Es mucho más rápido imshow				
    plt.yticks(plt.yticks()[0], [str(np.round(t*velUnit)) for t in plt.yticks()[0]]) 
    plt.ylabel("Velocity [km/s]",fontsize=fsize)
#    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
    plt.xticks(plt.xticks()[0], ["{:.2f}".format(t*estUnit) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]",fontsize=fsize)
    if(constantes[7] == JEANS):
#        plt.title("Jeans Instability $\\tau =$ {:d}".format(TAU),fontsize=fsize)
        plt.title("Phase Space Density $t =$ {:.2f} T".format(i*dt),fontsize=fsize) #Divido por 2 pues T = 2.
#        plt.ylim(constantes[0],constantes[1])
#        plt.xlim(constantes[2],constantes[3])
        #plt.clim(0,37e5) #Jeans
    elif(constantes[7] == GAUSS):
        #plt.title("Gaussian Initialization $\\tau =$ {:d}".format(TAU),fontsize=fsize)
#        plt.suptitle("$\\tau$ = {:d}".format(TAU),fontsize=fsize)
#        if(TAU == 0):
#            plt.suptitle("$\\tau$ = $\\infty$",fontsize=fsize)
        plt.title("Phase Space Density $t =$ {:.2f}".format(i*dt),fontsize=fsize)
        #plt.clim(0,30e5) #Gauss
        plt.xlim(constantes[2]/2,constantes[3]/2)
        plt.ylim(constantes[0]/1.5,constantes[1]/1.5)
    elif(constantes[7] == BULLET):
        #plt.title("Gaussian Initialization $\\tau =$ {:d}".format(TAU),fontsize=fsize)
#        plt.suptitle("$\\tau$ = {:d}".format(TAU),fontsize=fsize)
#        if(TAU == 0):
#            plt.suptitle("$\\tau$ = $\\infty$",fontsize=fsize)
        plt.ylim(constantes[2]/2,constantes[3]/2)
        plt.xlim(constantes[0],constantes[1])
        plt.title("Phase Space Density $t =$ {:.2f} ut".format(i*dt),fontsize=fsize)
        #plt.clim(0,27e5)

    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
    cbar.set_label("Mass density [$M_{\odot}$ / kpc  $\\frac{km}{s}$]",fontsize=fsize)
    
    plt.savefig("./images/phase{:d}.png".format(i), dpi = dpII)
    plt.clf()
    
    
    dens = np.loadtxt("./datFiles/density{:d}.dat".format(i))    
    plt.plot(x,dens)
    plt.xticks(plt.xticks()[0], ["{:.2f}".format(t*estUnit) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]",fontsize=fsize)
    plt.ylabel("Linear Density [$M_{\odot}$ / kpc]",fontsize=fsize)
    #plt.title("Density $\\tau =$ {:d}".format(TAU),fontsize=fsize)
    #plt.title("Density $t =$ {:.2f} ut".format(i*dt),fontsize=fsize)
    if(TAU == 0):
        plt.title("Density at {:.2f}, $\\tau \\rightarrow \\infty$".format(i*dt),fontsize=fsize)
    #plt.ylim(-0.75e9,6.85e10)#Gauss
#    plt.ylim(6e10,7e10)#Jeans
    plt.ylim(constantes[10]-0.05,constantes[10]+0.05)#Jeans
    #plt.xlim(constantes[0]/2,constantes[1]/2)
    plt.xlim(constantes[0],constantes[1])
    plt.savefig("./images/density{:d}.png".format(i), dpi = dpII)
    plt.clf()    
    
#    potential = np.loadtxt("./datFiles/potential{:d}.dat".format(i))
#    plt.plot(x,potential)
#    plt.ylabel("Potential [J /kg]",fontsize=fsize)
#    plt.xlim(-1.1,1.1)
#    plt.title("Potential $t =$ {:.2f} ut".format(i*dt),fontsize=fsize)
#    #plt.title("Potential $\\tau =$ {:d}".format(TAU),fontsize=fsize)
#    plt.ylim(-1.5e11,1.1e11)#Gauss
#    #plt.ylim(-1.6e11,1.1e11)#Jeans
#    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
#   # plt.yticks(plt.yticks()[0], [fmt(t*potUnit,1) for t in plt.yticks()[0]]) 
#    plt.xlabel("Position [kpc]",fontsize=fsize)
#    plt.savefig("./images/potential{:d}.png".format(i), dpi = dpII)
#    plt.clf()
#    
#    
#    acce = np.loadtxt("./datFiles/acce{:d}.dat".format(i))
#    plt.plot(x,acce)
#    plt.ylabel("Acceleration [kpc / Gy$^2$]",fontsize=fsize)
#    #plt.title("Acceleration $\\tau =-\\infty$",fontsize=fsize)
#    plt.ylim(-0.009,0.009)#Gauss
##    plt.ylim(-0.009,0.009)#Jeans
#    plt.xlim(-1.1,1.1) 
#    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
#    #plt.yticks(plt.yticks()[0], [fmt(t*acceUnit,1) for t in plt.yticks()[0]]) 
##    plt.ylim(np.min(acce)*1.1,np.max(acce)*1.1)
#    plt.title("Acceleration $t =$ {:.2f} ut".format(i*dt),fontsize=fsize)
#    plt.xlabel("Position [kpc]",fontsize=fsize)
#    
#    plt.savefig("./images/acce{:d}.png".format(i), dpi = dpII)
#    plt.clf()
    
f = open('plots', 'w+')
f.close()




