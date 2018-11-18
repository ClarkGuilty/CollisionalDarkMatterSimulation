#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 10:04:04 2018

@author: Javier Alejandro Acevedo Barroso
Script de Python para cálculos auxiliares.

"""
import numpy as np
#import pyqtgraph as pg
import matplotlib.pyplot as plt
import scipy as sc
import decimal as dec
D=dec.Decimal
    

#1.05457148 × 10-34 m2 kg / s
#retorna el valor de G en mis unidades
def unidades(x, mass, times):# x en megaparsecs, mass en masas solares y times en fracción de la edad del universo.
    x = D(x)
    mass = D(mass)
    times = D(times)
    x0 = D('3.0857e+22') #un megaparsec en metros.
    m0 = D('1.988e30')  #Masa solar en kg.
    t0 = D('13.772')*1000000000 #Edad del universo en años.
    t0 = t0*D('365.24')*24*60*60 #Ahora en segundos.
    G =  D('6.67408e-11')*np.power(x*x0,-3)*(m0*mass)*np.power(times*t0,2) #G en mis unidades.
    #hubble = 70/(x*x0*D('1e-3'))*(times*t0)*x #La constante de hubble en las unidades.
    #sv = 3e-26* np.power(x*x0*0.01,-3) * times*t0 
    h = D('1.05457148e-34')*np.power((x*x0),-2)/(m0*mass)*D(times*t0)
    #print("La constante de Hubble es: %f" % hubble)
    #print("La constante de Planck/2pi es:")
    #print(h)
    #print("La constante gravitacional es %f" % G)
    return G
    
#Retorna el valor de un metro / segundo en mis unidades.
def unidadesVel(x,mass,times):
    x = D(x)
    mass = D(mass)
    times = D(times)
    x0 = D('3.0857e+22' )
    t0 = D('13.772e9') * D('365.24')*24*60*60
    t = times*t0
    x = x*x0
    #sv = D('3e-26')*(x**(-3) ) *times*t0
    sv = x/ t
    return sv


#Retorna el valor de un metro / segundo^2 en mis unidades.
def unidadesAcce(x,mass,times):
    x = D(x)
    mass = D(mass)
    times = D(times)
    x0 = D('3.0857e+22' )
    t0 = D('13.772e9') * D('365.24')*24*60*60
    t = times*t0
    x = x*x0
    #sv = D('3e-26')*(x**(-3) ) *times*t0
    sv = x/ t**2
    return sv

#Retorna el valor de un metro^2 / segundo^2 en mis unidades.
def unidadesPot(x,mass,times):
    x = D(x)
    mass = D(mass)
    times = D(times)
    x0 = D('3.0857e+22' )
    t0 = D('13.772e9') * D('365.24')*24*60*60
    t = times*t0
    x = x*x0
    #sv = D('3e-26')*(x**(-3) ) *times*t0
    sv = x/ t
    return sv**2

#Calcula la masa de la materia oscura en mis unidades, massValue siendo la masa en eV.
def valorMasaDM(massValue, mass):
    mass = D(mass)
    massValue = D(massValue)
    val = D('1.783e-36')
    m0 = D('1.988e+30')
    return massValue*val/(m0*mass)
    
    
#Densidad de materia oscura
def valorDensidadMedia(x,mass,times):
    G = unidades(x,mass,times)
    
    x = D(x)
    mass = D(mass)
    times = D(times)
    x0 = D('3.0857e+22' )
    t0 = D('13.772e9') * D('365.24')*24*60*60
    hubble = D('67.4')/(x*x0*D('1e-3'))*(times*t0)*x #La constante de hubble en las unidades.
    t = times*t0
    x = x*x0
    #return 3*hubble**2/(8*D(np.pi)*G)*D('0.26')
    return 3*hubble**2/(8*D(np.pi)*G)*D('0.142')/(D(67.4)/100)**2
    
def valorSigmaV(x,mass,times):
    x = D(x)
    mass = D(mass)
    times = D(times)
    x0 = D('3.0857e+22' )
    t0 = D('13.772e9') * D('365.24')*24*60*60
    t = times*t0
    x = x*x0
    sv = D('3e-26')*np.power(x,-3) * t
    return sv

#massvalue en Electronvolts
def TAU(x,mass,times, massvalue):
    return valorMasaDM(massvalue,mass)/(valorDensidadMedia(x,mass,times)*
                       valorSigmaV(x,mass,times))

#x = 35e-3
#m = 0.1e12
#t= 4e-3

x = 50e-3
m = 0.1e12
t= 3e-3

#print("El valor de una unidad de velocidad es [km/s]}")
#print(unidadesVel(x,m,t)/1000)
#print("El valor de una unidad de potencial es [J/kg]")
#print(unidadesPot(x,m,t))
#print("El valor de una unidad de aceleración es [km/s²]")
#print(unidadesAcce(x,m,t)/1000)

newG = unidades(x,m,t) 
dmmass = 1000
#print("El 4*Pi*G = %f" % (4*np.pi*newG))
print("El G = %f" % (newG))
print("La densidad crítica %f" % (valorDensidadMedia(x,m,t)))
print("TAU es:")
print(TAU(x,m,t,700))

#print("La masa a usar es %f eV = %f" % (dmmass, unidadesMass(dmmass, 1e11)))
#print(unidadesMass(dmmass, 1e11))


def jeans(x,v,rho,sigma,A,k):
    return rho/np.sqrt(2*np.pi*sigma*sigma)*np.exp(-v*v/(2*sigma*sigma))*(1.0+A*np.cos(k*x))
    
#print(jeans(-0.57,-0.46, 10, 0.1, 4 , 2*np.pi))

