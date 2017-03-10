#!/usr/bin/env python
"""
Created on Thu Nov 24 22:29:44 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
cimport numpy as np
cimport cython
import scipy.constants as sc
from numpy import linalg as LA

@cython.boundscheck(False)
@cython.wraparound(False)

###############################################

#Defining Variables

cdef double AU, PC, R, Year, t_max, dt_max, e, n
cdef int Ng, Dump, t

#Size
AU = sc.astronomical_unit
PC = 206265*AU
R = 200*AU

#No.Of.groups
Ng = 1

#Dumping Number
Dump = 50

#Duration
Year = 365.26*(24*60*60)*(1.001)
t_max = (1*10**5)*Year
t = 0
dt_max = Year

#Initial Conditions

#Constants
e = 0.05*AU
n = 1*10**2

############################################

def M(int N):

    cdef int i
    cdef double M0
    cdef np.ndarray[np.double_t] Mass = np.zeros((N,))

    #Mass
    M0 = 1.989e30
    i = 0
    
    for i in range(N):
        M = M0    
        Mass[i] = np.array([M])
        
    return Mass

def GroupP(int Ng):

    cdef int i, O, Ns, S
    cdef double D, phi, X, Y, Z, C, As, R
    
    cdef np.ndarray[np.double_t, ndim=2] GroupPos = np.zeros((Ng,3))
    cdef np.ndarray[np.double_t] N = np.zeros((Ng,))
    
    C = 20000*AU
    As = 1000*AU
    R = 0.1*PC    
    
    i = -1
    O = -2
    Ns = 0
    
    while True:
        
        i = i + 1
        
        if i == (Ng):
            break
        
        elif O <= (PC-C)/C:
            
            S = np.random.randint(3, 6)
            
            N[i] = S
            
            Ns = Ns + S         
            
            O = O + 1    
            
            D = np.random.uniform(0, R**2)
            phi = np.random.uniform(0, 2*np.pi)           
            
            X = (O * C) + np.random.normal(C, As)
            Y = np.sqrt(D) * np.cos(phi)
            Z = np.sqrt(D) * np.sin(phi)
            
            GroupPos[i] = np.array([X,Y,Z])

        elif O > (PC-C)/C:
            
            S = np.random.randint(3, 6)
            
            N[i] = S
            
            Ns = Ns + S
            
            O = - 1
            
            D = np.random.uniform(0, R**2)
            phi = np.random.uniform(0, 2*np.pi)
            
            X = np.random.normal(0, As)
            Y = np.sqrt(D) * np.cos(phi)
            Z = np.sqrt(D) * np.sin(phi)
            
            GroupPos[i] = np.array([X,Y,Z])
            
    return GroupPos, Ns, N

def GroupV(int Ng):

    cdef int i
    cdef double Vx, Vy, Vz, C, Vgroup
    
    cdef np.ndarray[np.double_t, ndim=2] Vel = np.zeros((Ng,3))
    cdef np.ndarray[np.double_t] V = np.zeros((Ng,))

    i = 0

    for i in range(Ng):
        
        Vgroup = 0.2e3
        
        Vx, Vy, Vz = np.random.uniform(-1,1,3) 
        
        C = Vgroup / np.sqrt(Vx**2 + Vy**2 + Vz**2)
    
        Vel[i] = np.array([Vx * C, Vy * C, Vz * C])
        
        V[i] = LA.norm(Vel[i])
        
    return Vel
    
def PE(np.ndarray[np.double_t, ndim=2] Pos, np.ndarray[np.double_t] Mass, int e, int N):
      
    cdef int i, j    
    cdef double m, r
    
    cdef np.ndarray[np.double_t] Pe = np.zeros((N,))
    cdef double G = sc.gravitational_constant

    i = 0
    j = 0
    
    for i in range(0,N-1):
        for j in range(i+1,N):     
            
            r = Pos[i]-Pos[j]
            m = LA.norm(r)
            
            Pe[i] += -(G*Mass[i]*Mass[j])/(m+e) # Check Pe
            Pe[j] += -(G*Mass[j]*Mass[i])/(m+e)

    return Pe

def KE(np.ndarray[np.double_t, ndim=2] Vel, np.ndarray[np.double_t] Mass, int N):
    
    cdef double vi
    cdef int i
    cdef np.ndarray[np.double_t] Ke = np.zeros((N,))

    i = 0

    for i in range(0,N):
        vi = LA.norm(Vel[i])        
        Ke[i] = .5*Mass[i]*vi**2

    return Ke
  
def NormV( np.ndarray[np.double_t, ndim=2] Vel, np.ndarray[np.double_t, ndim=2] Pos, np.ndarray[np.double_t] Mass, int N):
    
    cdef double Ptot, Ktot, ERatio, l, V
    
    Ptot = np.sum(-PE(Pos, Mass, e, N))        
    
    Ktot = np.sum(KE(Vel, Mass, N))    
    
    ERatio = (2*Ktot)/Ptot
    
    l = np.random.uniform(0.9, 1)
    
    V = l*Vel/np.sqrt(ERatio)
    
    return V

########################################################

def IC(int Ns, int Ng, double R):

    V = []
    K = []
    P = []
    
    cdef int i, j, a, O
    cdef double X, Y, Z, MOD, Vx, Vy, Vz, GV
    
    cdef np.ndarray[np.double_t, ndim=2] Pos = np.zeros((Ns,3))
    cdef np.ndarray[np.double_t] Mass = np.zeros((Ns,))
    cdef np.ndarray[np.double_t, ndim=2] apos
    cdef np.ndarray[np.double_t, ndim=2] avel
    cdef np.ndarray[np.double_t] amass = np.zeros((N,))

    O = -1
    j = 0
    
    for j in range(Ng):
        i = -1
        a = int(N[j])
    
        amass = M(a)
        avel = np.zeros((a,3))
        apos = np.zeros((a,3))
        
        while i < a:
            i = i + 1
            O = O + 1
            
            X = np.random.uniform(-R, R)
            Y = np.random.uniform(-R, R)
            Z = np.random.uniform(-R, R)
            
            MOD = np.sqrt(X**2 + Y**2 + Z**2)
            
            if i == (N[j]):
                O = O - 1
                break
                 
            elif MOD <= R:
                apos[i,:] = np.array([X,Y,Z])
                Pos[O] = apos[i,:] + GroupPos[j]
                
                Vx, Vy, Vz = np.random.uniform(-1,1,3)
        
                avel[i] = np.array([Vx, Vy, Vz])
                
                Mass[O] = amass[i]
                
            elif MOD > R:
               i = i - 1
               O = O - 1
               
        GroupVel = GroupV(Ng)
    
        GV = NormV(avel, apos, amass, a) + GroupVel[j]
        
        K.append(np.sum(KE(GV, amass, a)))
        P.append(np.sum(PE(apos, amass, e, a)))
        
        V.append(GV)
    
    return Pos, V, Mass, K, P
        

def IV(V, int Ng, int Ns):
    
    cdef int j, k, c, O
    cdef np.ndarray[np.double_t, ndim=2] VF = np.zeros((Ns,3))
    
    O = -1
    j = 0
    
    for j in range(Ng):
        c = int(N[j])
        for k in range(c):
            O = O + 1
            
            VF[O] = V[j][k,:]
            
    return VF

#############################################
 
GroupPos, Ns, N = GroupP(Ng)

Pos, V, Mass, KinE, PotE = IC(Ns, Ng, R)

Vel = IV(V, Ng, Ns)
 
 
