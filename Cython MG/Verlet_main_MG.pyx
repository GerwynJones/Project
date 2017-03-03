#!/usr/bin/env python
"""
Created on Tue Oct 18 09:11:38 2016

@author: Admin
"""

from __future__ import division
cimport cython
import numpy as np 
import scipy.constants as sc
from numpy import linalg as LA

import pyximport
pyximport.install()

import Verlet_IC_MG

##################################################

def Verp(double oVel, double oPos, double dt, double acc):
    cdef double Pos
    "Position:"
    Pos = oPos + oVel*dt + .5*acc*dt**2
    return Pos
    
    
def Verv(double Pos, double Mass, double oVel, double dt, double acc, int e):
    cdef double anew, pe, Vel
    "Velocities:"
    anew, pe = Acc(double Pos, double Mass, double oVel, int e)
    Vel = oVel + .5*(acc + anew)*dt
    return Vel


def Acc(double Pos, double Mass, double Vel, int e):
    
    cdef int G, i, j
    cdef np.ndarray[np.double_t, ndim=2] acc
    cdef np.ndarray[np.double_t, ndim=1] Pe

    "Acceleration:"
    acc = np.zeros((Ns,3) dtype=np.double)
    Pe = np.zeros(Ns)
    G = sc.gravitational_constant
    
    for i in range(0,Ns-1):
        for j in range(i+1,Ns):
            
            r = Pos[i] - Pos[j]
            m = LA.norm(r)
            F = (G/(m+e)**3)*r # check (m+e) part
            
            acc[i] += -F*Mass[j]
            acc[j] += F*Mass[i]
            Pe[i] += -(G*Mass[i]*Mass[j])/(m+e) # Check PE
            Pe[j] += -(G*Mass[j]*Mass[i])/(m+e)

    return acc, Pe


def KE(double Vel, double Mass):
    
    cdef double vi
    cdef np.ndarray[np.double_t, ndim=2] Ke

    Ke = np.zeros(Ns)
    
    for i in range(0,Ns):
        vi = LA.norm(Vel[i])        
        Ke[i] = .5*Mass[i]*vi**2

    return Ke

############################################

def Loop(t, t_max) :
    
    P = []; A = []; E = []
    T = []; dT = []    
    
    cdef int O
    O = 0
    
    while t < t_max:   
        O = O + 1
        
        acc, Pe = Acc(Pos, Mass, Vel, e)
        
        Ke = KE(Vel,Mass)
        
        cdef double a, dt_grav, oPos, oVel, TE
        
        a = (LA.norm(acc, axis = 1))
        
        dt_grav =  np.min([dt_max, np.sqrt((2*n*e)/np.max(a))])
    
        "Verlet Method"
        
        oPos = Pos; oVel = Vel
    
        Pos = Verp(oVel, oPos, dt_grav, acc)
        Vel = Verv(Pos, Mass, oVel, dt_grav, acc, e)
    
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
        
        return P, A, T, dT, E