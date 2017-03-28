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
    """Position:"""
    Pos = oPos + (oVel * dt) + (0.5 * acc * dt ** 2)
    return Pos


def Verletv(Pos, Mass, oVel, dt, acc, e, Ns):
    """"Velocities:"""
    newacc, unused = Acc(Pos, Mass, e, Ns)
    Vel = oVel + (0.5 * (acc + newacc) * dt)
    return Vel


def Acc(Pos, Mass, e, Ns):
    """Acceleration:"""
    acc = np.zeros((Ns, 3))
    Pe = np.zeros(Ns)
    G = sc.gravitational_constant

    for i in xrange(0, Ns - 1):
        for j in xrange(i + 1, Ns):
            
            r = Pos[i] - Pos[j]
            magr = np.sqrt(inner1d(r, r))
            
            F = (G / ((magr + e) ** 3)) * r

            acc[i] += -F * Mass[j]
            acc[j] += F * Mass[i]

            Pe[i] += -(G * Mass[j] * Mass[i]) / (magr + e)

    return acc, Pe
    
    
def KE(Vel, Mass, Ns):
    """Kinetic Energy"""
    Ke = np.zeros(Ns)

    for i in xrange(0, Ns):
        
        vi = np.sqrt(inner1d(Vel[i], Vel[i]))
        Ke[i] = .5 * Mass[i] * vi ** 2

    return Ke

#############################################

P = []
E = []
T = []
dT = []
O = 0

############################################

while t < t_max:
    O = O + 1
    """ Calculating acceleration """
    acc, Pe = Acc(Pos, Mass, e, Ns)

    """ Energy Calculation """
    Ke = KE(Vel, Mass, Ns)
    TE = Ke + Pe    # Sum of star energies
    # Potential energy found with the acceleration function

    """ Calculating time-step """
    a = np.sqrt(inner1d(acc, acc))

    dt_grav = np.min([dt_max, np.sqrt((2 * eta * e) / np.max(a))])

    """ Verlet Method """
    oPos = Pos  # for previous steps with verlet
    oVel = Vel  # for previous steps with verlet

    Pos = Verletp(oVel, oPos, dt_grav, acc)
    Vel = Verletv(Pos, Mass, oVel, dt_grav, acc, e, Ns)

    """ Time """
    t = t + dt_grav

    if O == Dump:
        """ Dump Data into file """
        P.append(Pos)
        T.append(t + dt_grav)
        dT.append(dt_grav)
        E.append(TE)

        print ((t/t_max)*100)
        O = 0

    if t >= t_max:
        break

print t
