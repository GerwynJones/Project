#!/usr/bin/env python
"""
Created on Thu Nov 24 22:29:44 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
cimport numpy as np 
import scipy.constants as sc
from numpy import linalg as LA

@cython.boundscheck(False)
@cython.wraparound(False)

###############################################

" Defining Variables "

cdef double AU, PC, R, Year, t_max, dt_max, e, n
cdef int Ng, Dump, t

" Size "
AU = sc.astronomical_unit
PC = 206265*AU
R = 200*AU

" No.Of.groups "
Ng = 1

" Dumping Number"
Dump = 50

" Duration "
Year = 365.26*(24*60*60)*(1.001)
t_max = (1*10**5)*Year
t = 0
dt_max = Year

" Initial Conditions "

" Constants "
e = 0.05*AU
n = 1*10**2

############################################

def M(int N):

    cdef int i, M0
    cdef np.ndarray[np.double_t] Mass = np.zeros((N,), dtype=np.double)

    " Mass "
    M0 = 1.989e30
    
    for i in range(N):
        M = M0    
        Mass[i] = np.array([M])
        
    return Mass

def GroupP(int Ng):

    cdef int i, O, Ns
    cdef double S, D, phi, X, Y, Z, C, A, R
    
    cdef np.ndarray[np.double_t, ndim=2] GrouPos = np.zeros((Ng,3), dtype=np.double)   
    cdef np.ndarray[np.double_t] N = np.zeros((Ng,), dtype=np.double)
    
    C = 20000*AU
    A = 1000*AU    
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
            
            X = (O * C) + np.random.normal(C, A)
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
            
            X = np.random.normal(0, A)
            Y = np.sqrt(D) * np.cos(phi)
            Z = np.sqrt(D) * np.sin(phi)
            
            GroupPos[i] = np.array([X,Y,Z])
            
    return GroupPos, Ns, N

def GroupV(int Ng):

    cdef int i, Vgroup
    cdef double Vx, Vy, Vz, C
    
    cdef np.ndarray[np.double_t, ndim=2] Vel = np.zeros((Ng,3), dtype=np.double) 
    cdef np.ndarray[np.double_t] V = np.zeros((Ng,), dtype=np.double)
    
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
    
    cdef np.ndarray[np.double_t] Pe = np.zeros((N,), dtype=np.double)
    cdef double G = sc.gravitational_constant  
    
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
    cdef np.ndarray[np.double_t] Ke = np.zeros((N,), dtype=np.double)

    for i in range(0,N):
        vi = LA.norm(Vel[i])        
        Ke[i] = .5*Mass[i]*vi**2

    return Ke
  
def NormV( np.ndarray[np.double_t, ndim=2] Vel, np.ndarray[np.double_t, ndim=2] Pos, np.ndarray[np.double_t] Mass, int N):
    
    cdef double Ptot, Ktot, A, l, V
    
    Ptot = np.sum(-PE(Pos, Mass, e, N))        
    
    Ktot = np.sum(KE(Vel, Mass, N))    
    
    A = (2*Ktot)/Ptot
    
    l = np.random.uniform(0.9, 1)
    
    V = l*Vel/np.sqrt(A)    
    
    return V

########################################################

def IC(int Ns, int Ng, double R):

    V = []
    K = []
    P = []
    
    cdef int i, j, a, O
    cdef double mass, X, Y, Z, MOD, Vx, Vy, Vz, GroupVel, v
    
    cdef np.ndarray[np.double_t, ndim=2] Pos = np.zeros((Ns,3), dtype=np.double) 
    cdef np.ndarray[np.double_t] Mass = np.zeros((Ns,), dtype=np.double)    
    
    O = -1
    
    for j in range(Ng):
        i = -1
        a = int(N[j])
        
        cdef np.ndarray[np.double_t, ndim=2] pos = np.zeros((a,3), dtype=np.double) 
        cdef np.ndarray[np.double_t, ndim=2] vel = np.zeros((a,3), dtype=np.double)
    
        mass = M(a)
        
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
                pos[i] = np.array([X,Y,Z])
                Pos[O] = pos[i] + GroupPos[j]
                
                Vx, Vy, Vz = np.random.uniform(-1,1,3)
        
                vel[i] = np.array([Vx, Vy, Vz])
                
                Mass[O] = mass[i]
                
            elif MOD > R:
               i = i - 1
               O = O - 1
               
        GroupVel = GroupV(Ng)
    
        v = NormV(vel, pos, mass, a) + GroupVel[j]
        
        K.append(np.sum(KE(v, mass, a)))
        P.append(np.sum(PE(pos, mass, e, a))) 
        
        V.append(v)
    
    return Pos, V, Mass, K, P
        

def IV(V, int Ng, int Ns):
    
    cdef int j, k, a, O
    cdef double Vel
    cdef np.ndarray[np.double_t, ndim=2] Vel = np.zeros((Ns,3), dtype=np.double) 
    
    O = -1
    
    for j in range(Ng):
        a = int(N[j])
        for k in range(a):
            O = O + 1
            
            Vel[O] = V[j][k,:]
            
    return Vel

#############################################
 
GroupPos, Ns, N = GroupP(Ng)           

Pos, V, Mass, KinE, PotE = IC(Ns, Ng, R)

Vel = IV(V, Ng, Ns)
 
 
