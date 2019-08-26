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
L = constantes[1] - constantes[0]

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

#kj0 = np.fft.fftshift(np.loadtxt("./datFiles/powerSeries{:d}.dat".format(1)))[1049]
kj0=1

power = np.abs(np.loadtxt("./datFiles/powerSeries0.dat"))
freqs = np.fft.fftfreq(len(power), d = x[1]-x[0])


def findArg(array, value = 4*np.pi*2, epsilon = 0.1):
      heh = np.arange(len(power))
      while (len(heh[np.logical_and(array > value - epsilon, array < value + epsilon )]) == 0):
            epsilon += 0.1
      rta = heh[np.logical_and(freqs > value - epsilon, freqs < value + epsilon )]
      #print(rta[0])
      return rta[0]

fig, ax = plt.subplots()


#    ax.hlines(y=0.6, xmin=0.0, xmax=1.0, color='b')
power = np.abs(np.fft.fftshift(np.loadtxt("./datFiles/powerSeries0.dat")))
freqs = np.fft.fftshift(np.fft.fftfreq(len(power), d = x[1]- x[0])*L)

ax.set_yscale('log')
#ax.plot(freqs[power>1e-6],power[power>1e-6]**2)
ax.axhline(y=1.0, xmin=0.0, xmax=100.0, color='b', linewidth = 2)
ax.axvline(x=4*np.pi, color='r', linewidth = 0.5)
#    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
ax.set_xlabel("Position [kpc]",fontsize=fsize)
ax.set_ylabel("$P(k,t)/P(k_j,0)$",fontsize=fsize)
ax.set_title("Power spectrum $t =$ {:.2f} T".format(0*dt/2),fontsize=fsize)

fig.savefig("./images/powerSeries{:d}.png".format(0), dpi = dpII)
ax.cla()

kj = 2*np.pi/L
minylim = 1e-9
maxylim = 1e3
N = 1
nkj = findArg(freqs, value = kj)
#power = np.loadtxt("./datFiles/powerSeries1.dat")
power = np.fft.fftshift(np.loadtxt("./datFiles/powerSeries0.dat")) 
#power = np.abs(power[nkj-N:nkj+N])
#normalization = power.sum() /len(power[power!=0])
normalization = power[nkj]**2

for i in range(0,int(constantes[6])): 
      
    ax.cla()
    freqs = np.fft.fftshift(np.fft.fftfreq(len(power), d = x[1]- x[0])*L)
    power = np.abs(np.fft.fftshift(np.loadtxt("./datFiles/powerSeries{:d}.dat".format(i))))
    power[findArg(freqs,value = 0)] = 0
#    Ncentral = 1024
#    N = Ncentral-10
#    Nf = Ncentral+10
#    plt.plot(freqs[N:Nf],power[N:Nf])
#    ax.set_xscale('log')
    ax.set_yscale('log')
#    ax.plot(freqs[power>1e-6],power[power>1e-6]/power[1049])
    ax.axhline(y=1.0, xmin=0.0, xmax=10.0, color='b', linewidth = 1)
    ax.axvline(x=kj, color='r', linewidth = 0.5)
    
#    print(freqs.shape, (power/normalization).shape)
    ax.scatter(freqs[freqs != 0],power[freqs != 0]**2/normalization, s= 5)
#    ax.scatter(freqs[nkj-N:nkj+N],power[nkj-N:nkj+N]/normalization, s= 1)
#    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
    ax.set_xlabel("kL",fontsize=fsize)
    ax.set_ylabel("$|P(k,t)/P(k_j,0)|^2$",fontsize=fsize)
    #plt.title("Density $\\tau =$ {:d}".format(TAU),fontsize=fsize)
    ax.set_title("Power spectrum $t =$ {:.2f} T".format(i*dt/2),fontsize=fsize)
    #plt.ylim(-0.75e9,0.5e10)#Gauss
    #plt.ylim(-0.75e9,7e10)#Jeans
    ax.set_ylim(minylim, maxylim)
    ax.set_xlim(0, np.max(freqs))#
    fig.savefig("./images/powerSeries{:d}.png".format(i), dpi = dpII)

    






