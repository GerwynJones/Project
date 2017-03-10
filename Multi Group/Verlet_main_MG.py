#!/usr/bin/env python
"""
Created on Tue Oct 18 09:11:38 2016

@author: Admin
"""

from __future__ import division
import numpy as np
import scipy.constants as sc
from numpy import linalg as LA
import h5py

#from Verlet_IC_MG import *

##  File no.
Q = str( 1 )

# Position
with h5py.File('IC_No'+Q+'/Position_No'+Q+'.h5', 'r') as hf:
    Pos = hf['Position_Data'][:]

# Velocity
with h5py.File('IC_No'+Q+'/Velocity_No'+Q+'.h5', 'r') as hf:
    Vel = hf['Velocity_Data'][:]
    
# Mass
with h5py.File('IC_No'+Q+'/Mass_No'+Q+'.h5', 'r') as hf:
    Mass = hf['Mass_Data'][:]

# NGroup
with h5py.File('IC_No'+Q+'/Ng_No'+Q+'.h5', 'r') as hf:
    Ng = hf['Ng_Data'][:]

# Ns
with h5py.File('IC_No'+Q+'/Ns_No'+Q+'.h5', 'r') as hf:
    Ns = hf['Ns_Data'][:]
    
# N
with h5py.File('IC_No'+Q+'/N_No'+Q+'.h5', 'r') as hf:
    N = hf['N_Data'][:]
    
# Dump
with h5py.File('IC_No'+Q+'/Dump_No'+Q+'.h5', 'r') as hf:
    Dump = hf['Dump_Data'][:]
    
# T_max
with h5py.File('IC_No'+Q+'/Tmax_No'+Q+'.h5', 'r') as hf:
    t_max = hf['Tmax_Data'][:]
    
# T
with h5py.File('IC_No'+Q+'/t_No'+Q+'.h5', 'r') as hf:
    t = hf['t_Data'][:]
    
# dT
with h5py.File('IC_No'+Q+'/dt_No'+Q+'.h5', 'r') as hf:
    dt_max = hf['dt_Data'][:]
    
# eta
with h5py.File('IC_No'+Q+'/eta_No'+Q+'.h5', 'r') as hf:
    eta = hf['eta_Data'][:]
    
# e
with h5py.File('IC_No'+Q+'/e_No'+Q+'.h5', 'r') as hf:
    e = hf['e_Data'][:]
    

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
