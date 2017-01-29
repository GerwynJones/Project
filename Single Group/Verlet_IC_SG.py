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
R = (200/np.sqrt(3))*AU

" No.Of.Stars "

N = 3 

" Mass "
M0 = 1.989e30
M1 = 10*M0
M2 = 2*M0
M3 = 5*M0

mass = np.array([M1,M2,M3])

" Duration "

Year = 365.26*(24*60*60)*(1.001)
t_max = 1000*Year; t = 0; dt_max = Year/10

" Initial Conditions "

" Constants "

e = 0.05*AU; n = 1e-3; c = 1000

############################################

pos = np.zeros((N,3))
vel = np.zeros((N,3))

for i in range(N):
    
    X, Y, Z = np.random.uniform(-1,1,3)
    
    pos[i] = np.array([X*R,Y*R,Z*R])

    Vx, Vy, Vz = np.random.uniform(-1,1,3)

    vel[i] = np.array([Vx * c, Vy * c, Vz * c])

def PE(pos, mass, e):
    
    pe = np.zeros((N,1))
    G = sc.gravitational_constant    
    
    for i in range(0,N-1):
        for j in range(i+1,N):     
            
            r = pos[i]-pos[j]
            m = LA.norm(r)
            
            pe[i] += -(G*mass[i]*mass[j])/(m+e) # Check PE
            pe[j] += -(G*mass[j]*mass[i])/(m+e)

    return pe


def KE(vel, mass):
    
    ke = np.zeros((N,1))
    
    for i in range(0,N):
        vi = LA.norm(vel[i])        
        ke[i] = .5*mass[i]*vi**2

    return ke

Pos1,Pos2,Pos3 = pos
V1,V2,V3 = vel
P1,P2,P3 = -PE(pos, mass, e)
K1,K2,K3 = KE(vel, mass)

def normv(vel, pos, mass, PE):
    
    Ptot = np.sum(-PE(pos, mass, e))        
    
    Ktot = np.sum(KE(vel, mass))    
    
    A = Ktot/Ptot
    
    v = vel/np.sqrt(A)    
    
    return v, Ktot, Ptot
    
v, Ktot, Ptot = normv(vel, pos, mass, PE)    
    
   
    
    
    
    
    