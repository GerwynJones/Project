#!/usr/bin/env python
"""
Created on Thu Nov 24 22:59:34 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from Verlet_main_SG import *

plt.close('all')

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
plt.plot(A[:,0],A[:,1],A[:,2], label="1"); plt.plot(B[:,0],B[:,1],B[:,2], label="2") ; plt.plot(C[:,0],C[:,1],C[:,2], label="3"); plt.plot(D[:,0],D[:,1],D[:,2], label="4")
plt.legend(loc = 'best')

plt.figure(2)
plt.plot(T, Ea/Esum, label="1"), plt.plot(T, Eb/Esum, label="2"); plt.plot(T, ((Ea + Eb + Ec + Ed)/Esum), color = 'black') ; plt.plot(T, Ec/Esum, label="3"); plt.plot(T, Ed/Esum, label="4")
plt.legend(loc = 'best')
plt.savefig('Graph of Energy SG.png', bbox_inches='tight')

fig = plt.figure(3)
plt.plot(A[:,0],A[:,1], label="1"); plt.plot(B[:,0],B[:,1], label="2") ; plt.plot(C[:,0],C[:,1], label="3"); plt.plot(D[:,0],D[:,1], label="3")

plt.ylabel(r'Distance $(AU)$')
plt.xlabel(r'Distance $(AU)$')
plt.savefig('Graph of SG.png', bbox_inches='tight') 
plt.legend(loc = 'best')
#plt.xlim(-6e14,6e14); plt.ylim(-6e14,6e14)

plt.figure(4)
plt.plot(T, acceleration)

plt.show()

T_max = np.max(dT) 
T_min = np.min(dT)
T_Ratio = T_max/T_min

