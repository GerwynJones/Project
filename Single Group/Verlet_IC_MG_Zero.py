#!/usr/bin/env python
"""
Created on Thu Nov 24 22:29:44 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
import scipy.constants as sc
from numpy import linalg as LA

###############################################

" Defining Variables "

" Size "
AU = sc.astronomical_unit
PC = 206265*AU
R = 200*AU

" No.Of.Stars "
N = 11

" Mass "
M0 = 1.989e30
Mass = np.zeros(N)

for i in range(N):
    M = M0    
    Mass[i] = np.array([M])

" Duration "

Year = 365.26*(24*60*60)*(1.001)
t_max = 1e5*Year; t = 0; dt_max = Year

" Initial Conditions "

" Constants "

e = 0.05*AU; n = 1e-3

############################################

def GroupPos(N):
    
    Pos = np.zeros((N,3))
    
    C = 20000*AU
    
    A = 1000*AU    
    
    R = 0.1*PC    
    
    i = -1
    a = -2
    
    while True:
        
        i = i + 1
        
        if i == (N):
            break
        
        elif Pos[i-1,0] <= (PC-C):
            a = a + 1    
            
            D = np.random.uniform(0, R**2)
            phi = np.random.uniform(0, 2*np.pi)           
            
            X = (a * C) + np.random.normal(C, A)
            Y = np.sqrt(D) * np.cos(phi)
            Z = np.sqrt(D) * np.sin(phi)
            
            Pos[i] = np.array([X,Y,Z])

            
        elif Pos[i-1,0] > (PC-C):
            a = - 1
            
            D = np.random.uniform(0, R**2)
            phi = np.random.uniform(0, 2*np.pi)
            
            X = np.random.normal(0, A)
            Y = np.sqrt(D) * np.cos(phi)
            Z = np.sqrt(D) * np.sin(phi)
            
            Pos[i] = np.array([X,Y,Z])

            
    return Pos


def GroupVel(N):
    
    Vel = np.zeros((N,3))
    V = np.zeros(N)
    
    for i in range(N):
        
        Vgroup = 0.2e3
        
        Vx, Vy, Vz = np.random.uniform(-1,1,3) 
        
        C = Vgroup / np.sqrt(Vx**2 + Vy**2 + Vz**2)
    
        Vel[i] = np.array([Vx * C, Vy * C, Vz * C])
        
        V[i] = LA.norm(Vel[i])
        
    return Vel
    
    
#Vel = GroupVel(N)
Pos = GroupPos(N)


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.plot(Pos[:,0],Pos[:,1],Pos[:,2], marker = '.', linestyle = 'None')
plt.legend(loc = 'best')


plt.figure()
plt.plot(Pos[:,0], Pos[:,1], marker = '.', linestyle = 'None')
plt.grid()

plt.show()
