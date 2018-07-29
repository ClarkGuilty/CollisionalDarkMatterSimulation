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
newG = unidades(20e-3,1e12,3e-3) #Galaxia 
dmmass = 1000
print("El 4*Pi*G = %f" % (4*np.pi*newG))
#print("La masa a usar es %f eV = %f" % (dmmass, unidadesMass(dmmass, 1e11)))
print(unidadesMass(dmmass, 1e11))


def jeans(x,v,rho,sigma,A,k):
    return rho/np.sqrt(2*np.pi*sigma*sigma)*np.exp(-v*v/(2*sigma*sigma))*(1.0+A*np.cos(k*x))
    
print(jeans(-0.57,-0.46, 10, 0.1, 4 , 2*np.pi))
