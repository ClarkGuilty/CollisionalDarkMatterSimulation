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
    a, b = '{:.4e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

constantes = np.loadtxt("constants.dat", usecols = 1)
TAU = int(constantes[8])

#x = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
dpi = 500

densidadTheo = np.loadtxt('./datFiles/density0.dat').T
potTheo = np.loadtxt('./datFiles/potential0.dat').T
potReal = np.loadtxt('./datFiles/potential1.dat').T
#potReal = potReal - potReal[0,0]

plt.imshow(densidadTheo)
plt.title("Densidad cumple la ecuación de poisson del potTeorico")
cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
plt.savefig("densidadTeorica.png",dpi=dpi)

plt.figure()
plt.imshow(potTheo)
plt.title("potTeorico V = cos(0.5*pi*x)cos(0.5*pi*y)")
cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
plt.savefig("PotencialTeorico.png",dpi=dpi)

plt.figure()
plt.imshow(potReal)
plt.title("potCalculado a partir de la densidad")
cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
plt.savefig("potReal.png",dpi=dpi)


#accex = np.loadtxt('./datFiles/acce0x.dat').T
#accey = np.loadtxt('./datFiles/acce0y.dat').T
#
#gtheox = np.loadtxt('./datFiles/gtheox.dat').T
#gtheoy = np.loadtxt('./datFiles/gtheoy.dat').T

#lapl = (-4*potTheo[0:128,0:128]+potTheo[-1:127,1:129]+potTheo[1:129,1:129]+potTheo[-1:127,-1:127]+potTheo[-1:127,1:129])

def darX(inx):
    return -1.0+2.0/128*inx

def darY(iny):
    return -1.0+2.0/128*iny

dista = (int) (0.2*TAU)
mina = (int) (TAU//2 ) - dista 
maxa = (int) ( TAU//2+dista) 

#diff = (potTheo[mina:maxa,mina:maxa]-potReal[mina:maxa,mina:maxa])/potTheo[mina:maxa,mina:maxa]
diff = potTheo-potReal

plt.figure()
plt.imshow(diff)
plt.title("potTeorico/potCalculado")
cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
plt.savefig("potCalpotTeo.png", dpi=dpi)


realx = np.loadtxt('./datFiles/realx.dat').T
realy = np.loadtxt('./datFiles/realy.dat').T
calcx = np.loadtxt('./datFiles/calcx.dat').T
calcy = np.loadtxt('./datFiles/calcy.dat').T

fastShow(realx, title="realx")
fastShow(realy, title="realy")
fastShow(calcx, title="accex Calc - Teorica")
fastShow(calcy, title="accey Calc - Teorica")


#potT = potTheo.T-potTheo
#i,j = 3,3
#stencil = np.array([potTheo[i+1,j+1],potTheo[i-1,j+1],potTheo[i+1,j-1],potTheo[i-1,j-1]])
#stencil2 = np.array([potTheo[i,j-1],potTheo[i,j+1],potTheo[i+1,j],potTheo[i-1,j]])
#stencil3 = np.array([potTheo[i+1,j],potTheo[i,j+1]])
#print(i,j,potTheo[i+1,j+1],potTheo[i-1,j+1],potTheo[i+1,j-1],potTheo[i-1,j-1], stencil.sum()-4*potTheo[i,j])
#def rad(x,y):
#    return x*x+y*y
#x = np.linspace(-1,1,128)
#y = np.linspace(-1,1,128)
#
#values = np.zeros((128,128))
#values = rad(x,y)

#ratax = accex/gtheox
#ratay = accey/gtheoy
#
#plt.figure()
#plt.imshow(accey)
#plt.imshow(accex)
#plt.imshow(gtheoy)
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#plt.figure()
#plt.imshow(accey)
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#plt.figure()
#plt.imshow(ratax)
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#plt.figure()
#plt.imshow(np.log(ratay))
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))


#
    
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
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
#for i in range(128):
#    print(-1.0+i*2/127)



























