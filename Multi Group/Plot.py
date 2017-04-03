#!/usr/bin/env python
"""
Created on Thu Mar  2 23:50:10 2017

@author: Admin
"""
from __future__ import division
import numpy as np
import scipy.constants as sc
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import h5py

Q = str(3)

Ns = 5; Ng = 2; N = np.array([2, 3])

P = np.zeros(1000)

Energy = np.zeros((Ns, len(P)))
Esum = np.zeros((Ng, len(P)))
Time = np.zeros(len(P))

O = np.zeros(Ng + 1)
for k in xrange(int(Ns)):
    plt.figure()
    l = str(k+1)
    for i in xrange(N[k]):
        O[k] = O[k-1] + N[k]
        i = i + O[k-1]
        j = i + 1

        plt.plot(Time, Energy[i, :], label='star %s' % j)
    plt.ylabel("Energy (Ratio)")
    plt.xlabel("Time (yrs)")
    plt.legend(loc='best')
    plt.savefig('/home/gerwyn/Documents/Project-Large-Files/Graph.png', bbox_inches='tight')


