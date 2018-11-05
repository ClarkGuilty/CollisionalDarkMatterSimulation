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
import matplotlib
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 10})
plt.rcParams['image.cmap'] = 'plasma'
dpiT = 300
fsize = 16

constantes = np.loadtxt("./datFiles/constants.dat", usecols = 1)
#print(constantes)
tamano = int(constantes[14])
Nt = int(constantes[-1])
wMin = constantes[1]

figure = plt.figure(figsize=(7,5))

for i in range(len(constantes)):
    print(i,constantes[i])

def darV(i):
    return wMin+tamano*i

def fastShow(image, title="none",clim=None, clabel=None, saveN=None):
    plt.clf()
    plt.imshow(image,extent=[constantes[0],constantes[3],constantes[0],constantes[3]],norm=matplotlib.colors.LogNorm(vmin = clim[0], vmax = clim[1]), interpolation='nearest', aspect='auto')
#    plt.imshow(image,extent=[constantes[0],constantes[3],constantes[0],constantes[3]],interpolation='nearest', aspect='auto')
    cbar = plt.colorbar()
    if(clim != None):
        #plt.clim(clim[0],clim[1])
        print("heh")
    if(clabel != None):
        cbar.set_label(clabel,fontsize=fsize)
    plt.yticks(plt.yticks()[0], [str(np.round(t*50)) for t in plt.yticks()[0]])
    plt.ylabel("Position [kpc]",fontsize=fsize)
    plt.xticks(plt.xticks()[0], [str(np.round(t*50)) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]",fontsize=fsize)
    plt.title(title)
    plt.savefig("./images/"+saveN+".png",dpi=dpiT)




def fmt(x, pos):
    a, b = '{:.4e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)





densidadXY0  = np.loadtxt('./datFiles/densXY0.dat').T
densidadYZ0  = np.loadtxt('./datFiles/densYZ0.dat').T
densidadXZ0  = np.loadtxt('./datFiles/densXZ0.dat').T

#densidadXY1 = np.loadtxt('./datFiles/densXY1.dat').T
#densidadYZ1  = np.loadtxt('./datFiles/densYZ1.dat').T
#densidadXZ1  = np.loadtxt('./datFiles/densXZ1.dat').T


fastShow(densidadXY0, title="Corte Z = 0 de la densidad. Nt = {:d}".format(0),clim=[1e-4,5],clabel="Mass density [$M_{\odot}$ / kp$c^3$]",saveN="XY0" )
#fastShow(densidadYZ0, "YZ0")
#fastShow(densidadXZ0, "XZ0")

#fastShow(densidadXY1, "XY1")
#fastShow(densidadYZ1, "YZ1")
#fastShow(densidadXZ1, "XZ1")
#densidadXY = np.loadtxt("./datFiles/densXY{:d}.dat".format(i)).T
#plt.figure()
imagenes = tamano//8
#for i in range(0,tamano,tamano//imagenes):
#    potXY0 = np.loadtxt("./datFiles/pot0XY{:d}.dat".format(i)).T
    #potXY1 = np.loadtxt("./datFiles/pot1XY{:d}.dat".format(i)).T
 #   accexXY0 = np.loadtxt("./datFiles/accex0XY{:d}.dat".format(i)).T
    #accexXY1 = np.loadtxt("./datFiles/accex1XY{:d}.dat".format(i)).T

#    fastShow( potXY0 , "potXY{:d}".format(i) )
    #fastShow( accexXY0 , "accex0XY{:d}".format(i) )
   # fastShow( accexXY1 , "accex1XY{:d}".format(i) )
    

for i in range(1,Nt):
	densidadXY = np.loadtxt('./datFiles/densXY{:d}.dat'.format(i)).T
	densidad2 = np.log(densidadXY)
#	densidad2 = densidad2[densidad2 != -inf
	#densidadYZ = np.loadtxt('./datFiles/densYZ{:d}.dat'.format(i)).T
	#densidadXZ = np.loadtxt('./datFiles/densXZ{:d}.dat'.format(i)).T
	fastShow(densidadXY,title="Corte Z = 0 de la densidad. Nt = {:d}".format(i), clim=[1e-4,5],clabel="Mass density [$M_{\odot}$ / kp$c^3$]",saveN="XY{:d}".format(i) )
	#fastShow(densidadYZ, "YZ{:d}".format(i))
	#fastShow(densidadXZ, "XZ{:d}".format(i))






















