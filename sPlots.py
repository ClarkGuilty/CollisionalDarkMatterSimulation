# -*- coding: utf-8 -*-
"""


@author: Javier Alejandro Acevedo Barroso
Script de Python para comparar versión colisional con no colisional.

"""
import numpy as np
#import pyqtgraph as pg
import matplotlib.pyplot as plt
import scipy as sc


constantes = np.loadtxt("constants.dat", usecols = 1)
x = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
#inF = np.loadtxt("inF.dat")
#outF = np.loadtxt("outF0.dat")
#outF1 = np.loadtxt("outF1.dat")
#oI = np.loadtxt("oI.dat")
#oR = np.loadtxt("oR.dat")
#acce = np.loadtxt("acce.dat")


#        
#plt.pcolormesh(dat)



for i in range(int(constantes[6])):
    col = np.loadtxt("./col/grid{:d}.dat".format(i)).T
    nocol = np.loadtxt("./nocol/grid{:d}.dat".format(i)).T
    diff = nocol - col
    

    plt.imshow(diff, extent=[constantes[0],constantes[1],constantes[2],constantes[3]]) #Es mucho más rápido imshow
    plt.title("Phase-Space Difference")
    plt.colorbar()
#    plt.xticks(plt.yticks()[0], [str(t) for t in plt.yticks()[0]])
    plt.xticks(plt.xticks()[0], [str(t/constantes[1]) for t in plt.xticks()[0]])
#    plt.yticks(plt.yticks()[0], [str(t/constantes[3]) for t in plt.yticks()[0]])
    plt.xticks(plt.xticks()[0], [str(t) for t in plt.xticks()[0]])
    plt.savefig("./dif/phase{:d}.png".format(i))
    plt.clf()

    dcol = np.loadtxt("./col/density{:d}.dat".format(i))
    dnocol = np.loadtxt("./nocol/density{:d}.dat".format(i))
    ddiff = dnocol - dcol
    plt.plot(x,ddiff)
    plt.plot((0, 0), (-1, 1), 'k-')
    plt.title("Density Difference")
    plt.savefig("./dif/density{:d}.png".format(i))
    plt.clf()

    pcol = np.loadtxt("./col/potential{:d}.dat".format(i))
    pnocol = np.loadtxt("./nocol/potential{:d}.dat".format(i))
    pdiff = pnocol - pcol
    plt.plot(x,pdiff)
    plt.title("Potential Difference")
    plt.savefig("./dif/potential{:d}.png".format(i))
    plt.clf()

    acol = np.loadtxt("./col/acce{:d}.dat".format(i))
    anocol = np.loadtxt("./nocol/acce{:d}.dat".format(i))
    adiff = anocol - acol
    plt.plot(x,adiff)
    plt.title("Acceleration Difference")
    plt.savefig("./dif/acce{:d}.png".format(i))
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
