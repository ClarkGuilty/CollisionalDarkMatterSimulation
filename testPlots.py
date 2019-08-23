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

#kj0 = np.fft.fftshift(np.loadtxt("./datFiles/powerSeries{:d}.dat".format(1)))[1049]
kj0=1

power = np.abs(np.loadtxt("./datFiles/powerSeries0.dat"))
freqs = np.fft.fftfreq(len(power), d = 1.0/1024)
heh = np.arange(len(power))
heh[freqs == 4*np.pi*2]
heh[np.logical_and(freqs > 4*np.pi*2 - 0.3, freqs < 4*np.pi*2 + 0.3 )]
fig, ax = plt.subplots()


#    ax.hlines(y=0.6, xmin=0.0, xmax=1.0, color='b')
power = np.abs(np.fft.fftshift(np.loadtxt("./datFiles/powerSeries0.dat")))
freqs = np.fft.fftshift(np.fft.fftfreq(len(power), d = 1.0/1024)*2)
#    Ncentral = 1024
#    N = Ncentral-10
#    Nf = Ncentral+10
#    plt.plot(freqs[N:Nf],power[N:Nf])
#    ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(freqs[power>1e-6],power[power>1e-6])
ax.axhline(y=1.0, xmin=0.0, xmax=100.0, color='b', linewidth = 2)
ax.axvline(x=4*np.pi, color='r', linewidth = 0.5)
#    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
ax.set_xlabel("Position [kpc]",fontsize=fsize)
ax.set_ylabel("Linear Density [$M_{\odot}$ / kpc]",fontsize=fsize)
#plt.title("Density $\\tau =$ {:d}".format(TAU),fontsize=fsize)
ax.set_title("Power spectrum $t =$ {:.2f} T".format(0*dt/2),fontsize=fsize)
#plt.ylim(-0.75e9,0.5e10)#Gauss
#plt.ylim(-0.75e9,7e10)#Jeans
#    ax.set_ylim(1e-7, 1e3)
#ax.set_xlim(1e-1, 1e3)
fig.savefig("./images/powerSeries{:d}.png".format(0), dpi = dpII)
ax.cla()

minylim = 1e-7

N = 1
nkj = 1037
power = np.fft.fftshift(np.loadtxt("./datFiles/powerSeries1.dat")) 
#power = np.abs(power[nkj-N:nkj+N])
#normalization = power.sum() /len(power[power!=0])

for i in range(1,int(constantes[6])): 
    
#    ax.axhline(y=0.5, xmin=0.0, xmax=1.0, color='r')
#    ax.hlines(y=0.6, xmin=0.0, xmax=1.0, color='b')
    freqs = np.fft.fftshift(np.fft.fftfreq(2048, d = 1.0/1024)*2)
    power = np.abs(np.fft.fftshift(np.loadtxt("./datFiles/powerSeries{:d}.dat".format(i))))
    
#    Ncentral = 1024
#    N = Ncentral-10
#    Nf = Ncentral+10
#    plt.plot(freqs[N:Nf],power[N:Nf])
    ax.set_xscale('log')
    ax.set_yscale('log')
#    ax.plot(freqs[power>1e-6],power[power>1e-6]/power[1049])
    ax.axhline(y=1.0, xmin=0.0, xmax=10.0, color='b', linewidth = 1)
    ax.axvline(x=4*np.pi, color='r', linewidth = 0.5)
    
#    print(freqs.shape, (power/normalization).shape)
    ax.scatter(freqs,power, s= 1)
#    ax.scatter(freqs[nkj-N:nkj+N],power[nkj-N:nkj+N]/normalization, s= 1)
#    plt.xticks(plt.xticks()[0], [str(t*estUnit) for t in plt.xticks()[0]])
    ax.set_xlabel("Position [kpc]",fontsize=fsize)
    ax.set_ylabel("Linear Density [$M_{\odot}$ / kpc]",fontsize=fsize)
    #plt.title("Density $\\tau =$ {:d}".format(TAU),fontsize=fsize)
    ax.set_title("Power spectrum $t =$ {:.2f} T".format(i*dt/2),fontsize=fsize)
    #plt.ylim(-0.75e9,0.5e10)#Gauss
    #plt.ylim(-0.75e9,7e10)#Jeans
    ax.set_ylim(minylim, 1e3)
#    ax.set_xlim(1e-1, 1e3)
    fig.savefig("./images/powerSeries{:d}.png".format(i), dpi = dpII)
    ax.cla()
    


f = open('plots', 'w+')
f.close()




