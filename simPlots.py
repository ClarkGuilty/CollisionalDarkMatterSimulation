# -*- coding: utf-8 -*-
"""


@author: Javier Alejandro Acevedo Barroso
Script de Python para la visualización de la simulación.

"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sc

#dat = np.loadtxt("grid.dat")
density = np.loadtxt("density.dat")
constantes = np.loadtxt("constants.dat", usecols = 1)
inF = np.loadtxt("inF.dat")
outF = np.loadtxt("outF0.dat")
outF1 = np.loadtxt("outF1.dat")
oI = np.loadtxt("oI.dat")
oR = np.loadtxt("oR.dat")


  
#        
#plt.pcolormesh(dat)
#plt.savefig("phase.png")
#plt.figure()
#xf = np.linspace(0,1-1/constantes[4],int(constantes[4])) //Espacio de frecuencias senoidal
#x = np.linspace(0,2*np.pi, int(constantes[4]))           //Espacio x senoidal.
#plt.plot(x,density)
#plt.savefig("densidad.png")
#plt.figure()
#sns.distplot(density, kde=False, rug=True);






x = np.linspace(constantes[0], constantes[1], int(constantes[4]))
#plt.plot(x,inF)
#plt.plot(x,oR, color = 'black')
#plt.savefig("potencial.png")
#plt.scatter(x,density, color ='g')

#Verifica que el arreglo final luego de un loop sea igual al inicial.
#diffPorc = sum(np.abs(oR -density)/density)*100
#print diffPorc
#print "La diferencia porcentual entre el arreglo original y el arreglo un loop de fourier después es: "+ str(np.floor(diffPorc*10000)/10000)+"%"



#plt.plot(x,np.sin(x))
#plt.plot(x,outF1)
plt.plot(density)
plt.savefig("density.png")
h = plt.figure()
plt.plot(outF)
plt.savefig("Fourier.png")
h = plt.figure()
plt.plot(x,oR)
plt.savefig("Potencial.png")
#plt.plot(x,outM)
diferencia = oR-density

#simps(simps(z, y), x)