# -*- coding: utf-8 -*-
"""


@author: Javier Alejandro Acevedo Barroso
Script de Python para comparar versión colisional con no colisional.

"""
import numpy as np
#import pyqtgraph as pg
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib.ticker as ticker
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

constantes = np.loadtxt("constants.dat", usecols = 1)
x = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
TAU = int(constantes[8])
#inF = np.loadtxt("inF.dat")
#outF = np.loadtxt("outF0.dat")
#outF1 = np.loadtxt("outF1.dat")
#oI = np.loadtxt("oI.dat")
#oR = np.loadtxt("oR.dat")
#acce = np.loadtxt("acce.dat")

def fmt(x, pos):
    #a = '{:2}'.format(x)
    a = "{0:.0f}".format(x)
    return r'${:2} \% $'.format(a)
#        
#plt.pcolormesh(dat)


figu = plt.gcf()
dpII = 150
for i in range(int(constantes[6])):
    col = np.loadtxt("./col/grid{:d}.dat".format(i)).T
    nocol = np.loadtxt("./nocol/grid{:d}.dat".format(i)).T
    col = col/sum(sum(col))
    nocol = nocol/sum(sum(nocol))
    
    dif = nocol-col
    dif = dif*100
    
#    for aq in range(2048):
#        for aw in range(2048):
#            if(nocol[aq][aw]>0):
#                dif[aq][aw] = (nocol[aq][aw]-col[aq][aw])/nocol[aq][aw]*100
#            else:
#                dif[aq][aw] = 0
    
    
#    for z in range(2047):
#        for y in range(2047):
#            if(nocol[z][y] != 0):
#                diff[z][y] = (1.0)/nocol[z][y]
#            else:
#                diff[z][y] = 0


    plt.imshow(dif, extent=[constantes[0],constantes[1],constantes[2],constantes[3]]) #Es mucho más rápido imshow
    #plt.yticks(plt.yticks()[0], [str(np.round(t*473,2)) for t in plt.yticks()[0]]) 
    plt.ylabel("Velocity [km/s]")
    plt.xticks(plt.xticks()[0], [str(t*200) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]")
    plt.title("Gauss Comparison $(\\tau = 0$) - ($\\tau =$ {:d})".format(TAU))
        
    #cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
    #cbar = plt.colorbar(format = '%.0f%%')
    cbar=plt.colorbar(cmap='plasma')
    
    #cbar.set_ticks([0, .2, .4, .6, .8, 1])
    #cbar.set_ticklabels(['0', '20%','40%','60%', '80%', '100%'])
    cbar.set_label("Percentual difference")
    cbar.set_clim(-0.0004,0.0005)
    plt.savefig("./dif/phase{:d}.png".format(i), dpi=dpII)
    plt.clf()



#------------------------------------

    dcol = np.loadtxt("./col/density{:d}.dat".format(i))
    dnocol = np.loadtxt("./nocol/density{:d}.dat".format(i))

    
    dcol = dcol/sum(dcol)
    dnocol = dnocol/sum(dnocol)
    ddif = (dnocol - dcol)*100
#    for aq in range(2047):
#        if(dnocol[aq]>0):
#            ddif[aq] = (dnocol[aq]-dcol[aq]/dnocol[aq])*100
#        else:
#            ddif[aq]= 0
    plt.plot(x,ddif)
    #plt.plot((0, 0), (-1, 1), 'k-')
    plt.title("Density Comparison $(\\tau = 0$) - $\\tau =$ {:d}".format(TAU))
    plt.xticks(plt.xticks()[0], [str(t*200) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]")
    
    plt.ylabel("Porcentual difference")
    plt.ylim(-0.14, 0.20)
    plt.savefig("./dif/density{:d}.png".format(i), dpi=dpII)
    plt.clf()
    
    
#----------------------------------
    
    
    pcol = np.loadtxt("./col/potential{:d}.dat".format(i))
    pnocol = np.loadtxt("./nocol/potential{:d}.dat".format(i))

    pcol = pcol/sum(pcol)
    pnocol = pnocol/sum(pnocol)
    pdif = (pnocol - pcol)*100
    
    
    plt.plot(x,pdif)
    plt.xticks(plt.xticks()[0], [str(t*200) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]")
    plt.ylabel("Porcentual difference")
    plt.ylim(-0.0006, 0.001)
    plt.title("Potential Comparison $(\\tau = 0$) - $\\tau =$ {:d}".format(TAU))
    plt.savefig("./dif/potential{:d}.png".format(i), dpi=dpII)
    plt.clf()


#-------------------------------------------------
    
    
    acol = np.loadtxt("./col/acce{:d}.dat".format(i))
    anocol = np.loadtxt("./nocol/acce{:d}.dat".format(i))
    acol = acol/sum(acol)
    anocol = anocol/sum(anocol)
    adif = (anocol - acol)*100
    plt.plot(x,adif)
    plt.xticks(plt.xticks()[0], [str(t*200) for t in plt.xticks()[0]])
    plt.xlabel("Position [kpc]")
    plt.ylabel("Porcentual difference")
    plt.title("Acceleration Difference $(\\tau = 0$) - $\\tau =$ {:d}".format(TAU))
    plt.savefig("./dif/acce{:d}.png".format(i), dpi=dpII)
    plt.clf()

#    potential = np.loadtxt("./datFiles/potential{:d}.dat".format(i))
#    plt.plot(x,potential)
#    plt.savefig("./images/potential{:d}.png".format(i))
#    plt.clf()
#    acce = np.loadtxt("./datFiles/acce{:d}.dat".format(i))
#    plt.plot(x,acce)
#    plt.savefig("./images/acce{:d}.png".format(i))
#    plt.clf()
    
#xf = np.linspace(0,1-1/constantes[4],int(constantes[4])) #Espacio de frecuencias senoidal
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
