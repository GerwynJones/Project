#!/usr/bin/env python
"""
Created on Tue Oct 18 09:11:38 2016

@author: Admin
"""

from __future__ import division
import numpy as np
import scipy.constants as sc
from numpy import linalg as LA

#  File no.
Q = str( 1 )

npzfile = np.load('IC/IC_No'+Q+'.npz')

Pos = npzfile['arr_0']; Vel = npzfile['arr_1']; Mass = npzfile['arr_2']
Ns = npzfile['arr_3']; Ng = npzfile['arr_4']; N = npzfile['arr_5']
e = npzfile['arr_6']; eta = npzfile['arr_7']
PC = npzfile['arr_8']; R = npzfile['arr_9']
t_max = npzfile['arr_10']; dt_max = npzfile['arr_11']; t = npzfile['arr_12']
Dump = npzfile['arr_13']

##################################################

def Verp(oVel, oPos, dt, acc):
    "Position:"
    Pos = oPos + (oVel * dt) + (0.5 * acc * dt ** 2)
    return Pos


def Verv(Pos, Mass, oVel, dt, acc, e, Ns):
    "Velocities:"
    anew, pe = Acc(Pos, Mass, e, Ns)
    Vel = oVel + (0.5 * (acc + anew) * dt)
    return Vel


def Acc(Pos, Mass, e, Ns):
    "Acceleration:"
    acc = np.zeros((Ns, 3))
    Pe = np.zeros(Ns)
    G = sc.gravitational_constant

    for i in range(0, Ns - 1):
        for j in range(i + 1, Ns):
            r = Pos[i] - Pos[j]
            m = LA.norm(r)
            F = (G / (m + e) ** 3) * r  # check (m+e) part

            acc[i] += -F * Mass[j]
            acc[j] += F * Mass[i]
            Pe[i] += -(G * Mass[i] * Mass[j]) / (m + e)  # Check PE
            Pe[j] += -(G * Mass[j] * Mass[i]) / (m + e)

    return acc, Pe


def KE(Vel, Mass, Ns):
    "Kinetic Energy"
    Ke = np.zeros(Ns)

    for i in range(0, Ns):
        vi = LA.norm(Vel[i])
        Ke[i] = .5 * Mass[i] * vi ** 2

    return Ke

#############################################

P = []
A = []
E = []
T = []
dT = []
O = 0

############################################

while t < t_max:
    O = O + 1

    acc, Pe = Acc(Pos, Mass, e, Ns)

    Ke = KE(Vel, Mass, Ns)

    a = LA.norm(acc, axis=1)

    dt_grav = np.min([dt_max, np.sqrt((2 * eta * e) / np.max(a))])

    "Verlet Method"

    oPos = Pos
    oVel = Vel

    Pos = Verp(oVel, oPos, dt_grav, acc)
    Vel = Verv(Pos, Mass, oVel, dt_grav, acc, e, Ns)

    t = t + dt_grav
    TE = Ke + Pe

    if O == Dump:
        """Dump Data into file"""

        P.append(Pos)
        A.append(a)
        T.append(t + dt_grav)
        dT.append(dt_grav)
        E.append(TE)

        print ((t/t_max)*100)
        O = 0

    if t >= t_max:
        break
