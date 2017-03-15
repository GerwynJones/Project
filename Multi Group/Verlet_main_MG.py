#!/usr/bin/env python
"""
Created on Tue Oct 18 09:11:38 2016

@author: Admin
"""

from __future__ import division
import numpy as np
import scipy.constants as sc
from numpy.core.umath_tests import inner1d


from Verlet_IC_MG import *

##################################################

def Verletp(oVel, oPos, dt, acc):
    "Position:"
    Pos = oPos + (oVel * dt) + (0.5 * acc * dt ** 2)
    return Pos


def Verletv(Pos, Mass, oVel, dt, acc, e, Ns):
    "Velocities:"
    anew = Acc(Pos, Mass, e, Ns)
    Vel = oVel + (0.5 * (acc + anew) * dt)
    return Vel


def Acc(Pos, Mass, e, Ns):
    "Acceleration:"
    acc = np.zeros((Ns, 3))
    G = sc.gravitational_constant

    for i in xrange(0, Ns - 1):
        for j in xrange(i + 1, Ns):
            
            r = Pos[i] - Pos[j]

            magr = np.sqrt(inner1d(r, r))
            
            F = (G / ((magr + e) ** 3)) * r

            acc[i] += -F * Mass[j]
            acc[j] += F * Mass[i]

    return acc

def PE(Pos, Mass, e, Ns):
    "Potential Energy"

    Pe = np.zeros(Ns)
    G = sc.gravitational_constant

    for i in xrange(0, Ns - 1):
        for j in xrange(i + 1, Ns):
            
            r = Pos[i] - Pos[j]

            magr = np.sqrt(inner1d(r, r))

            Pe[i] += -(G * Mass[j] * Mass[i]) / (magr + e)
            Pe[j] += -(G * Mass[i] * Mass[j]) / (magr + e)  # Check PE
            
    return Pe
    
    
def KE(Vel, Mass, Ns):
    "Kinetic Energy"
    Ke = np.zeros(Ns)

    for i in xrange(0, Ns):
        
        vi = np.sqrt(inner1d(Vel[i], Vel[i]))
        
        Ke[i] = .5 * Mass[i] * vi ** 2

    return Ke

#############################################

P = []
N1 = []; N2 = []
E = []
T = []
dT = []
O = 0

############################################

while t < t_max:
    " Calculating acceleration and time-step "
    
    O = O + 1

    acc = Acc(Pos, Mass, e, Ns)

    a = np.sqrt(inner1d(acc, acc))

    dt_grav = np.min([dt_max, np.sqrt((2 * eta * e) / np.max(a))])

    "Verlet Method"

    oPos = Pos
    oVel = Vel

    Pos = Verletp(oVel, oPos, dt_grav, acc)
    Vel = Verletv(Pos, Mass, oVel, dt_grav, acc, e, Ns)
    
    "Energy Calculation"
    
    Ke = KE(Vel, Mass, Ns)
    Pe = PE(Pos, Mass, e, Ns)
    
    TE = Ke + Pe    # SUM

    " Time "
    
    t = t + dt_grav

    if O == Dump:
        """Dump Data into file"""

        P.append(Pos)
        T.append(t + dt_grav)
        dT.append(dt_grav)
        E.append(TE)
        N1.append(Ke)
        N2.append(Pe)

        print ((t/t_max)*100)
        O = 0

    if t >= t_max:

        break

print t
