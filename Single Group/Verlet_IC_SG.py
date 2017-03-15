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
R = 200*AU

" No.Of.Stars "
N = 2

" Mass "
M0 = 1.989e30
Mass = np.zeros(N)

for i in range(N):
    M = M0    
    Mass[i] = np.array([M])

" Duration "

Year = 365.26*(24*60*60)*(1.001)
t_max = 1000*Year; t = 0; dt_max = Year/10

" Initial Conditions "

" Constants "

e = 0.05*AU; n = 1e-3

############################################

Pos = np.zeros((N,3))
Vel = np.zeros((N,3))

i = -1

while True:
    
    i = i + 1

    X = np.random.uniform(-R, R)
    Y = np.random.uniform(-R, R)
    Z = np.random.uniform(-R, R)
    
    if i == (N):
        break
    
    elif np.sqrt(X**2 + Y**2 + Z**2) <= R:

        Pos[i] = np.array([X,Y,Z])
        
        Vx, Vy, Vz = np.random.uniform(-1,1,3)

        Vel[i] = np.array([Vx, Vy, Vz])
        
    elif np.sqrt(X**2 + Y**2 + Z**2) > R:
       i = i - 1

#############################################
            
def PE(Pos, Mass, e):
    
    Pe = np.zeros((N,1))
    G = sc.gravitational_constant    
    
    for i in range(0,N-1):
        for j in range(i+1,N):     
            
            r = Pos[i]-Pos[j]
            m = LA.norm(r)
            
            Pe[i] += -(G*Mass[i]*Mass[j])/(m+e) # Check Pe
            Pe[j] += -(G*Mass[j]*Mass[i])/(m+e)

    return Pe

def KE(Vel, Mass):
    
    Ke = np.zeros((N,1))
    
    for i in range(0,N):
        vi = LA.norm(Vel[i])        
        Ke[i] = .5*Mass[i]*vi**2

    return Ke

def NormV(Vel, Pos, Mass, PE):
    
    Ptot = np.sum(-PE(Pos, Mass, e))        
    
    Ktot = np.sum(KE(Vel, Mass))    
    
    A = (2*Ktot)/Ptot
    
    V = Vel/np.sqrt(A)    
    
    return V, Ptot
    
Vel, Ptot = NormV(Vel, Pos, Mass, PE)
Ktot = np.sum(KE(Vel, Mass))    
    

   
    
    
    
    
    