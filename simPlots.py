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


JEANS = -137
GAUSS = -127

dat = np.loadtxt("./datFiles/grid0.dat").T
#density = np.loadtxt("density.dat")
constantes = np.loadtxt("constants.dat", usecols = 1)
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
dpII = 150
#figu.tight_layout()
#plt.pcolormesh(dat)
for i in range(int(constantes[6])):
    dat = np.loadtxt("./datFiles/grid{:d}.dat".format(i)).T
    plt.imshow(dat, extent=[constantes[0],constantes[1],constantes[2],constantes[3]], interpolation='none', aspect='auto') #Es mucho más rápido imshow
    #plt.pcolormesh(dat)
    #plt.pcolor(dat) Nunca usar para grandes grillas	
    plt.yticks(plt.yticks()[0], [str(np.round(t*473,2)) for t in plt.yticks()[0]]) 
    plt.ylabel("Velocity [km/s]")
    plt.xticks(plt.xticks()[0], [str(t*200) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]")
    if(constantes[7] == JEANS):
        plt.title("Jeans Instability $\\tau =$ {:d}".format(TAU))
    elif(constantes[7] == GAUSS):
        plt.title("Gaussian Initialization $\\tau =$ {:d}".format(TAU))
    cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
    cbar.set_label("Mass density $M_{\odot}$ / kpc  $\\frac{km}{s}$")
    plt.savefig("./images/phase{:d}.png".format(i))
    cbar.set_clim(0,1e5)
    
    cbar.set_label("Mass density [$M_{\odot}$ / kpc  $\\frac{km}{s}$]")
    plt.savefig("./images/phase{:d}.png".format(i), dpi = dpII)
    plt.clf()
    
    
    dens = np.loadtxt("./datFiles/density{:d}.dat".format(i))    
    #plt.xticks(str(plt.xticks()[0]*200), [str(t*200) for t in plt.xticks()[0]])
    plt.plot(x,dens)
    plt.xticks(plt.xticks()[0], [str(t*200) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]")
    plt.ylabel("Linear Density $M_{\odot}$/kpc")
    plt.title("Density")
    plt.ylim(-0.75e9,1.6e10)
    #plt.plot((0, 0), (-1, 1), 'k-')
    plt.savefig("./images/density{:d}.png".format(i), dpi = dpII)
    plt.clf()
    
    
    potential = np.loadtxt("./datFiles/potential{:d}.dat".format(i))
    plt.plot(x,potential)
    plt.ylabel("Potential [J /kg]")
    plt.title("Potential")
    plt.ylim(-6.6e10,-5.8e10)
    plt.xticks(plt.xticks()[0], [str(t*200) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]")
    plt.savefig("./images/potential{:d}.png".format(i), dpi = dpII)
    plt.clf()
    
    
    acce = np.loadtxt("./datFiles/acce{:d}.dat".format(i))
    plt.plot(x,acce)
    plt.ylabel("Acceleration [kpc / $(mYears)^2$]")
    plt.title("Acceleration")
    #plt.yticks(plt.yticks()[0], [str(t*2754463327) for t in plt.yticks()[0]])
    plt.xticks(plt.xticks()[0], [str(t*200) for t in plt.xticks()[0]])
    
    plt.ylim(-0.0007,0.0007)
    
    plt.xlabel(" [kpc]")
    plt.savefig("./images/acce{:d}.png".format(i), dpi = dpII)
    plt.clf()
    
xf = np.linspace(0,1-1/constantes[4],int(constantes[4])) #Espacio de frecuencias senoidal
#plt.plot(x,density)
#plt.savefig("densidad.png")
#plt.figure()
#sns.distplot(density, kde=False, rug=True);

#plt.plot(x,inF)
#plt.plot(x,oR, color = 'black')
#plt.savefig("potencial.png")
#plt.scatter(x,density, color ='g')

#Verifica que el arreglo final luego de un loop sea igual al inicial.
#diffPorc = sum(np.abs(oR -density)/density)*100
#print diffPorc
#print "La diferencia porcentual entre el arreglo original y el arreglo un loop de fourier después es: "+ str(np.floor(diffPorc*10000)/10000)+"%"


#h = plt.figure()
#plt.plot(outF)
#plt.savefig("Fourier.png")
#h = plt.figure()
#plt.plot(x,oR)
#plt.savefig("Potencial.png")
#h = plt.figure()
#plt.plot(x,acce)
#plt.savefig("Aceleracion.png")

#simps(simps(z, y), x)
