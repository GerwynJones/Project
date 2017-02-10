#!/usr/bin/env python
"""
Created on Thu Nov 24 22:59:34 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from Verlet_main_MG import *

plt.close('all')

Position = np.zeros((Ns,len(P),3))
Acceleration = np.zeros((Ns,len(P)))
Time = np.zeros(len(P))

for j in range(len(P)):
    for i in range(Ns):
        Position[i][j,:] = P[j][i,:] 
        Acceleration[i,j] = A[j][i]
        Time[j] = T[j]
        
        
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')          
        
for i in range(Ns):    
    plt.plot(Position[i][:,0],Position[i][:,1],Position[i][:,2]) 

plt.xlabel("x")
plt.legend(loc = 'best')


#plt.figure(2)
#plt.plot(T, Ea/Esum, label="1"), plt.plot(T, Eb/Esum, label="2"); plt.plot(T, ((Ea + Eb + Ec + Ed)/Esum), color = 'black') ; plt.plot(T, Ec/Esum, label="3"); plt.plot(T, Ed/Esum, label="4")
#plt.legend(loc = 'best')
#plt.savefig('Graph of Energy SG.png', bbox_inches='tight')


plt.figure(3)

for i in range(Ns):    
    plt.plot(Position[i][:,0],Position[i][:,1]) 

plt.ylabel(r'Distance $(AU)$')
plt.xlabel(r'Distance $(AU)$')
plt.savefig('Graph of 2D MG.png', bbox_inches='tight') 
plt.legend(loc = 'best')

plt.figure(2)

for i in range(Ns):    
    plt.plot(T, Acceleration[i,:])

plt.show()

T_max = np.max(dT) 
T_min = np.min(dT)
T_Ratio = T_max/T_min

