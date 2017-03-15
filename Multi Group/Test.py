#!/usr/bin/env python
"""
Created on Thu Mar  2 23:50:10 2017

@author: Admin
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.close('all')


def M(N, mean = np.log10(0.08), sigma = 0.69):

    """ Mass """
    M0 = 1.989e30
    Mass = np.zeros(N)

    C = np.log(10)/(0.158*np.sqrt(2*np.pi)*sigma)
    
    for i in range(N):
        
        M = np.random.lognormal(mean, sigma)*C
    
        Mass[i] = np.array([M*M0])

    return Mass
    

N = 1000
M0 = 1.989e30

Mass = M(N)

plt.figure()

plt.hist(Mass/M0, bins=np.logspace(-1, 1.5, 100))

plt.xscale('log')
plt.xlabel(r'$Log_{10} M$')
plt.ylabel(r'N')
plt.title("Histogram of mass of stars")


Mpdf  = np.linspace(0, 30, 400)    
pdf = 0.158*(np.exp(-(np.log10(Mpdf) - np.log10(0.08))**2. /(2.*(0.69**2.))) / (np.log(10.)*Mpdf))

plt.figure()

plt.loglog(Mpdf, pdf, linewidth=2, color='r')

plt.axis('tight')

plt.show()
