# -*- coding: utf-8 -*-
"""


@author: Javier Alejandro Acevedo Barroso
Script de Python para la visualización de la simulación.

"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sc

dat = np.loadtxt("grid.dat")
density = np.loadtxt("density.dat")
constantes = np.loadtxt("constants.dat", usecols = 1)


#for i in range(100):
#    for j in range(100):
#        dat[i][j] -= i*100+j
        
        
plt.pcolormesh(dat)
plt.savefig("phase.png")
plt.figure()
x = np.linspace(-1,1,int(constantes[4]))
plt.plot(x,density)
plt.savefig("densidad.png")
plt.figure()
sns.distplot(density, kde=False, rug=True);
#print dat.sum()/np.sqrt(np.pi)


#simps(simps(z, y), x)