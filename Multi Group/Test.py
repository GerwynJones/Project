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

def IMF(alpha, C, M):
    B = (C/(-alpha + 1))*(M**(-alpha + 1))
    return B

def M(N):

    """ Mass """
    M0 = 1.989e30
    Mass = np.zeros(N)

    alpha = 2.3

    M_Min = 1
    M_Max = 100

    N_Min = 1

    C = (N_Min*(-alpha + 1))/(M_Max**(-alpha + 1))

    N_Max = IMF(alpha, C, M_Min)

    O = -1

    while True:
        O = O + 1
        dN = np.random.uniform(0,N_Max)
        M = np.random.uniform(M_Min,M_Max)

        B = IMF(alpha, C, M)

        if O == N:
            break

        elif dN <= B:
            Mass[O]= np.array([M*M0])

        else:
            O = O - 1

    return Mass, C


N = 100000

Mass, C = M(N)

x = np.linspace(1,100, N)
y = IMF(2.3, C, x) + 1600


plt.figure()
plt.hist(Mass/(1.989e30)) # , bins=np.logspace(0, 2, 100))
plt.plot(x, y)

#plt.xscale('log')
plt.xlabel(r'$Log_{10} M$')
plt.ylabel(r'N')
plt.title("Histogram of mass of stars")


plt.show()
