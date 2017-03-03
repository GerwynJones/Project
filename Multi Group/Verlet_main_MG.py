#!/usr/bin/env python
"""
Created on Tue Oct 18 09:11:38 2016

@author: Admin
"""

from __future__ import division
import numpy as np 
import scipy.constants as sc
from numpy import linalg as LA

from Verlet_IC_MG import *

##################################################

def Verp(oVel, oPos, dt, acc):
    "Position:"
    Pos = oPos + oVel*dt + .5*acc*dt**2
    return Pos
    
    
def Verv(Pos, Mass, oVel, dt, acc, e, Ns):
     "Velocities:"
     anew, pe = Acc(Pos, Mass, oVel, e, Ns)
     Vel = oVel + .5*(acc + anew)*dt
     return Vel


def Acc(Pos, Mass, Vel, e, Ns):
    "Acceleration:"
    acc = np.zeros((Ns,3))
    Pe = np.zeros(Ns)
    G = sc.gravitational_constant
    
    for i in range(0, Ns-1):
        for j in range(i+1, Ns):
            
            r = Pos[i] - Pos[j]
            m = LA.norm(r)
            F = (G/(m+e)**3)*r # check (m+e) part
            
            acc[i] += -F*Mass[j]
            acc[j] += F*Mass[i]
            Pe[i] += -(G*Mass[i]*Mass[j])/(m+e) # Check PE
            Pe[j] += -(G*Mass[j]*Mass[i])/(m+e)

    return acc, Pe


def KE(Vel, Mass, Ns):
    
    Ke = np.zeros(Ns)
    
    for i in range(0,Ns):
        vi = LA.norm(Vel[i])        
        Ke[i] = .5*Mass[i]*vi**2

    return Ke

#############################################

P = []; A = []; E = []
T = []; dT = []

O = 0 

############################################

while t < t_max:   
    O = O + 1
    
    acc, Pe = Acc(Pos, Mass, Vel, e, Ns)
    
    Ke = KE(Vel, Mass, Ns)
    
    a = (LA.norm(acc, axis = 1))
    
    dt_grav =  np.min([dt_max, np.sqrt((2*n*e)/np.max(a))])

    "Verlet Method"
    
    oPos = Pos; oVel = Vel

    Pos = Verp(oVel, oPos, dt_grav, acc)
    Vel = Verv(Pos, Mass, oVel, dt_grav, acc, e, Ns)

    t += dt_grav
    TE = Ke + Pe

    if O == Dump :   
        """Dump Data into file"""
    
        P.append(Pos)
        A.append(a)
        T.append(t + dt_grav)
        dT.append(dt_grav)    
        E.append(TE)  
        
        O = 0
    
    if t == t_max:
        break