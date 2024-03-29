#!/usr/bin/env python
"""
Created on Thu Nov 24 22:29:44 2016

@author: Admin
"""
from __future__ import division
import numpy as np
import scipy.constants as sc
from numpy import linalg as LA
# import h5py

from Kroupa_IMF import *

###############################################

# Defining Variables
# Size
AU = sc.astronomical_unit
G = sc.gravitational_constant
PC = 206265*AU
R = 200*AU

# No.Of.groups
Ng = 10

# Dumping Number
Dump = 60

# Duration
Year = sc.Julian_year
t_max = 1e6*Year; t = 0; dt_max = 2*Year

# Initial Conditions
# Constants
e = 0.01*AU; eta = 25

# Type
Virial = 1/2
Cold = Virial/2
Hot = Virial*2

############################################

def M(N, Method):

    alpha = np.array([1.35, 2.35])
    
    M_min = 0.1
    M_max = 100

    """ Mass """
    M0 = 1.989e30

    M = Method(N, alpha, M_min, M_max)
    
    Mass = M*M0

    return Mass

def GroupP(Ng):
    
    GroupPos = np.zeros((Ng, 3))
    
    N = np.zeros(Ng)
    
    C = 20000*AU
    
    A = 1000*AU    
    
    R = 0.1*PC    
    
    i = -1
    O = -2
    Ns = 0
    
    while True:
        
        i = i + 1
        
        if i == Ng:
            break
        
        elif O <= (PC-2*C)/C:
            
            S = np.random.randint(3, 6)
            
            N[i] = S
            
            Ns = Ns + S         
            
            O = O + 1    
            
            D = np.random.uniform(0, R**2)
            phi = np.random.uniform(0, 2*np.pi)           
            
            X = (O * C) + np.random.normal(C, A)
            Y = np.sqrt(D) * np.cos(phi)
            Z = np.sqrt(D) * np.sin(phi)
            
            GroupPos[i] = np.array([X, Y, Z])

        elif O > (PC-2*C)/C:

            S = np.random.randint(3, 6)
            
            N[i] = S
            
            Ns = Ns + S
            
            O = - 1
            
            D = np.random.uniform(0, R**2)
            phi = np.random.uniform(0, 2*np.pi)
            
            X = np.random.normal(0, A)
            Y = np.sqrt(D) * np.cos(phi)
            Z = np.sqrt(D) * np.sin(phi)
            
            GroupPos[i] = np.array([X, Y, Z])
            
    return GroupPos, Ns, N

def GroupV(Ng):
    
    Vel = np.zeros((Ng, 3))
    V = np.zeros(Ng)
    
    for i in xrange(Ng):
        
        Vgroup = 0.2e3
        
        Vx, Vy, Vz = np.random.uniform(-1, 1, 3)
        
        C = Vgroup / np.sqrt(Vx**2 + Vy**2 + Vz**2)
    
        Vel[i] = np.array([Vx * C, Vy * C, Vz * C])
        
        V[i] = LA.norm(Vel[i])
        
    return Vel
    
def PE(Pos, Mass, e, N):
    
    Pe = np.zeros(N)

    for i in xrange(0, N - 1):
        for j in xrange(i + 1, N):     
            
            r = Pos[i]-Pos[j]
            m = LA.norm(r)
            
            Pe[i] += -(G*Mass[i]*Mass[j])/(m+e)  # Check Pe

    return Pe

def KE(Vel, Mass, N):
    
    Ke = np.zeros(N)
    
    for i in xrange(0, N):
        modv = LA.norm(Vel[i])
        Ke[i] = 0.5*Mass[i]*modv**2

    return Ke
  
def NormV(Vel, Pos, Mass, N, Type):
    
    Ptot = np.sum(-PE(Pos, Mass, e, N))        
    
    " This is not Virial as it doesnt take into account the applied group velocity "
    
    Ktot = np.sum(KE(Vel, Mass, N))
    
    Tot = (Ptot*Type)/Ktot
    
    V = Vel*np.sqrt(Tot)
    
    return V

########################################################

def IC(Ns, Ng, N, R, GroupPos, Method, Type):
    """ Creating the initial conditions with random points within a cylinder """    
    
    Pos = np.zeros((Ns, 3))
    Mass = np.zeros(Ns)
    V = []
    K = []
    P = []
    
    O = -1
    
    for j in xrange(Ng):
        i = -1
        a = int(N[j])
        
        apos = np.zeros((a, 3))
        avel = np.zeros((a, 3))

        """ Creating masses from either Salpeter or Kroupa"""
        amass = M(a, Method)
        
        while i < a:
            i = i + 1
            O = O + 1
            
            X = np.random.uniform(-R, R)
            Y = np.random.uniform(-R, R)
            Z = np.random.uniform(-R, R)
            
            if i == (N[j]):
                O = O - 1
                break
                 
            elif np.sqrt(X**2 + Y**2 + Z**2) <= R:
                apos[i] = np.array([X, Y, Z])
                Pos[O] = apos[i] + GroupPos[j]
                
                Vx, Vy, Vz = np.random.uniform(-1, 1, 3)
        
                avel[i] = np.array([Vx, Vy, Vz])
                
                Mass[O] = amass[i]
                
            elif np.sqrt(X**2 + Y**2 + Z**2) > R:
               i = i - 1
               O = O - 1
               
        GroupVel = GroupV(Ng)
    
        GV = NormV(avel, apos, amass, a, Type) + GroupVel[j]
        
        Velf = NormV(GV, apos, amass, a, Type)
        
        K.append(np.sum(KE(Velf, amass, a)))
        P.append(np.sum(PE(apos, amass, e, a)))
        
        V.append(Velf)
    
    return Pos, V, Mass, K, P, Type
        

def IV(V, Ng, Ns, N):
    """ Converting the Velocity array into something the loop can use """
    
    Vel = np.zeros((Ns, 3))
    
    O = -1
    
    for j in xrange(Ng):
        a = int(N[j])
        for k in xrange(a):
            O = O + 1
            
            Vel[O] = V[j][k, :]
            
    return Vel

#############################################
 
GroupPos, Ns, N = GroupP(Ng)

Pos, V, Mass, KinE, PotE, Type = IC(Ns, Ng, N, R, GroupPos, Kroupa, Virial)

Vel = IV(V, Ng, Ns, N)



# Earth-Sun IC

# Ns = 2; N = np.array([2])
#
# Mass = np.array([1.989e30, 5.972e24])
#
# Pos = np.zeros((Ns,3))
#
# Pos[1] = np.array([0,AU,0])
#
# Vel = np.zeros((Ns,3))
#
# V = 29754.7
#
# Vel[1] = np.array([V,0,0])
