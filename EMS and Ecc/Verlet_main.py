#!/usr/bin/env python
"""
Created on Tue Oct 18 09:11:38 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
import scipy.constants as sc
from numpy.core.umath_tests import inner1d

from Verlet_IC_EMS import *

#####################

def Verp(ovel, opos, dt, a):
    "Position:"
    pos = opos + ovel*dt + .5*a*dt**2
    return pos
    
def Verv(pos, mass, ovel, dt, a, e):
     "Velocities:"
     adt = acc(pos, mass, e)
     vel = ovel + .5*(a + adt)*dt
     return vel


def acc(pos, mass, e):
    "Acceleration:"
    a = np.zeros((N, 3))
    G = sc.gravitational_constant

    for i in xrange(0, N - 1):
        for j in xrange(i + 1, N):
            
            r = pos[i] - pos[j]

            magr = np.sqrt(inner1d(r, r))
            
            F = (G / ((magr + e) ** 3)) * r

            a[i] += -F * mass[j]
            a[j] += F * mass[i]

    return a

def PE(pos, mass, e):
    "Potential Energy"

    pe = np.zeros(N)
    G = sc.gravitational_constant

    for i in xrange(0, N - 1):
        for j in xrange(i + 1, N):
            
            r = pos[i] - pos[j]

            magr = np.sqrt(inner1d(r, r))

            pe[i] += -(G * mass[j] * mass[i]) / (magr + e)
#            pe[j] += -(G * mass[i] * mass[j]) / (magr + e)  # Check PE
            
    return pe


def KE(vel, mass):
    
    ke = np.zeros(N)
    
    for i in xrange(0, N):
        vi = np.sqrt(inner1d(vel[i],vel[i]))      
        ke[i] = .5*mass[i]*vi**2

    return ke

while t < t_max:   
    " Calculating acceleration and time-step "
    
    ac = acc(pos, mass, e)
    
    a_o = np.sqrt(inner1d(ac, ac))
    
    dt_grav = np.min([dt_max, np.sqrt((2*n*e)/np.max(a_o))])

    print(t/t_max)*100 
    T.append(t + dt_grav)

    
    "Verlet Method"
    opos = pos; ovel = vel

    pos = Verp(ovel, opos, dt_grav, ac)
    vel = Verv(pos, mass, ovel, dt_grav, ac, e)

    t += dt_grav
    
    ke = KE(vel,mass)
    pe = PE(pos, mass, e)
    
    """Dump pos into file"""
    a.append(pos[0])
    A = np.asarray(a)
    b.append(pos[1])
    B = np.asarray(b)
    c.append(pos[2])
    C = np.asarray(c)

    ea.append(pe[0] + ke[0])
    Ea = np.asarray(ea)
    eb.append(pe[1] + ke[1])
    Eb = np.asarray(eb)
    ec.append(pe[2] + ke[2])
    Ec = np.asarray(ec)

    Tsum.append(np.sum(ke + pe)) 
    Esum = np.asarray(Tsum)
    
    if t == t_max:
        break
    
