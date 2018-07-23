# -*- coding: utf-8 -*-
"""


@author: Javier Alejandro Acevedo Barroso
Script de Python para cálculos auxiliares.

"""
import numpy as np
#import pyqtgraph as pg
import matplotlib.pyplot as plt
import scipy as sc
import decimal as dec

D=dec.Decimal
    

def unidades(x, mass, times):# x en megaparsecs, mass en masas solares y times en fracción de la edad del universo.
    x0 = 3.0857e+22 #un megaparsec en metros.
    m0 = 1.988e30  #Masa solar en kg.
    t0 = 13.772*1000000000 #Edad del universo en años.
    t0 = t0*365.24*24*60*60 #Ahora en segundos.
    G =  6.67408e-11*np.power(x*x0,-3)*(m0*mass)*np.power(times*t0,2) #G en mis unidades.
    hubble = 70/(x*x0*1e-3)*(times*t0)*x #La constante de hubble en las unidades.
    sv = 3e-26* np.power(x*x0*0.01,-3) * times*t0 
    print("La constante de Hubble es: %f" % hubble)
    print("La constante gravitacional es %f" % G)
    return G
    

def unidadesDec(x,mass,times):
    x = D(x)
    mass = D(mass)
    times = D(times)
    x0 = D('3.0857e+20' )
    t0 = D('13.772e9') * D('365.24')*24*60*60
    x = x*x0
    sv = D('3e-26')*(x**(-3) ) *times*t0
    return sv

#Calcula la masa de la materia oscura en mis unidades, massValue siendo la masa en eV.
def unidadesMass(massValue, mass):
    mass = D(mass)
    massValue = D(massValue)
    val = D('1.783e-36')
    m0 = D('1.988e+30')
    return massValue*val/(m0*mass)
    
    
print(unidadesDec(20e-3,1e11,3e-3))

#newG = unidades(5,1e15,2e-1) #Clusters
newG = unidades(20e-3,1e11,3e-3) #Galaxia 
dmmass = 1000
print("El 4*Pi*G = %f" % (4*np.pi*newG))
#print("La masa a usar es %f eV = %f" % (dmmass, unidadesMass(dmmass, 1e11)))
print(unidadesMass(dmmass, 1e11))
#constantes = np.loadtxt("constants.dat", usecols = 1)
#x = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
#inF = np.loadtxt("inF.dat")
#outF = np.loadtxt("outF0.dat")
#outF1 = np.loadtxt("outF1.dat")
#oI = np.loadtxt("oI.dat")
#oR = np.loadtxt("oR.dat")
#acce = np.loadtxt("acce.dat")


#        
#plt.pcolormesh(dat)



#for i in range(int(constantes[6])):
#    dat = np.loadtxt("./datFiles/grid{:d}.dat".format(i))
#    plt.imshow(dat, extent=[constantes[0],constantes[1],constantes[2],constantes[3]]) #Es mucho más rápido imshow
#    #plt.pcolormesh(dat)
#    #plt.pcolor(dat) Nunca usar para grandes grillas	
#    plt.savefig("./images/phase{:d}.png".format(i))
#    plt.clf()
#dens = np.loadtxt("density.dat")
#plt.plot(x,dens)
#    plt.savefig("./images/density{:d}.png".format(i))
#    plt.clf()
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
