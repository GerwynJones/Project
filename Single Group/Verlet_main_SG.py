#!/usr/bin/env python
"""
Created on Tue Oct 18 09:11:38 2016

@author: Admin
"""

from __future__ import division
import numpy as np 
import scipy.constants as sc
from numpy import linalg as LA

from Verlet_IC_MG_Zero import *

##################################################

def Verp(oVel, oPos, dt, acc):
    "Position:"
    Pos = oPos + oVel*dt + .5*acc*dt**2
    return Pos
    
    
def Verv(Pos, Mass, oVel, dt, acc, e):
     "Velocities:"
     anew, pe = Acc(Pos, Mass, oVel, e)
     Vel = oVel + .5*(acc + anew)*dt
     return Vel


def Acc(Pos, Mass, Vel, e):
    
    acc = np.zeros((N,3))
    Pe = np.zeros((N,1))
    G = sc.gravitational_constant
    
    for i in range(0,N-1):
        for j in range(i+1,N):
            
            r = Pos[i] - Pos[j]
            m = LA.norm(r)
            F = (G/(m+e)**3)*r # check (m+e) part
            
            acc[i] += -F*Mass[j]
            acc[j] += F*Mass[i]
            Pe[i] += -(G*Mass[i]*Mass[j])/(m+e) # Check PE
            Pe[j] += -(G*Mass[j]*Mass[i])/(m+e)

    return acc, Pe


def KE(Vel, Mass):
    
    Ke = np.zeros((N,1))
    
    for i in range(0,N):
        vi = LA.norm(Vel[i])        
        Ke[i] = .5*Mass[i]*vi**2

    return Ke

#############################################

a = []; Ta = []
b = []; Tb = []
c = []; Tc = []
d = []; Td = []

Tsum = []

T = []; dT = []


t = 0; ac = []

############################################

while t < t_max:   
    
    acc, Pe = Acc(Pos, Mass, Vel, e)
    
    Ke = KE(Vel,Mass)
    
    norm_a = (LA.norm(acc, axis = 1)); ac.append(norm_a); acceleration = np.asarray(ac)
    
    dt_grav =  np.min([dt_max, np.sqrt((2*n*e)/np.max(norm_a))])

    print(t/t_max)*100 
    T.append(t + dt_grav)
    dT.append(dt_grav)    
    
    "Verlet Method"
    oPos = Pos; oVel = Vel

    Pos = Verp(oVel, oPos, dt_grav, acc)
    Vel = Verv(Pos, Mass, oVel, dt_grav, acc, e)

    t += dt_grav
    
    """Dump Pos into file"""
    a.append(Pos[0]/AU)
    A = np.asarray(a)
    b.append(Pos[1]/AU)
    B = np.asarray(b)
    c.append(Pos[2]/AU)
    C = np.asarray(c)
    d.append(Pos[3]/AU)
    D = np.asarray(d)
    """Dump energies into file"""
    Ta.append(Ke[0]/2 + Pe[0])
    Ea = np.asarray(Ta)
    Tb.append(Ke[1]/2 + Pe[1])
    Eb = np.asarray(Tb)
    Tc.append(Ke[2]/2 + Pe[2])
    Ec =  np.asarray(Tc)
    Td.append(Ke[3]/2 + Pe[3])
    Ed =  np.asarray(Td)
    Tsum.append(Pe[0] + Ke[0]/2 + Pe[1] + Ke[1]/2 + Pe[2] + Ke[2]/2 + Pe[3] + Ke[3]/2)
    Esum =  np.asarray(Tsum)
    
    if t == t_max:
        break
    
