#!/usr/bin/env python
"""
Created on Thu Nov 24 22:29:44 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
import scipy.constants as sc
from numpy import linalg as LA


###############################################

# Defining Variables

# Size
AU = sc.astronomical_unit
PC = 206265*AU
R = 200*AU

# No.Of.groups
Ng = 10

# Dumping Number
Dump = 50

# Duration
Year = 365.26*(24*60*60)
t_max = 1e6*Year; t = 0; dt_max = 5*Year

# Initial Conditions

# Constants
e = 0.1*AU; eta = 5

############################################

def IMF(alpha, C, M):
    B = (C/(-alpha + 1))*(M**(-alpha + 1))
    return B

def M(N):

    """ Mass """
    M0 = 1.989e30
    Mass = np.zeros(N)

    alpha = 2.3

    M_Min = 1
    M_Max = 100

    N_Min = 1

    C = (N_Min*(-alpha + 1))/(M_Max**(-alpha + 1))

    N_Max = IMF(alpha, C, M_Min)

    O = -1

    while True:
        O = O + 1

        dN = np.random.uniform(0,N_Max)
        M = np.random.uniform(M_Min,M_Max)

        B = IMF(alpha, C, M)

        if O == N:
            break

        elif dN <= B:
            Mass[O]= np.array([M*M0])

        else:
            O = O - 1

    return Mass

def GroupP(Ng):
    
    GroupPos = np.zeros((Ng,3))
    
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
            
            GroupPos[i] = np.array([X,Y,Z])

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
            
            GroupPos[i] = np.array([X,Y,Z])
            
    return GroupPos, Ns, N

def GroupV(Ng):
    
    Vel = np.zeros((Ng,3))
    V = np.zeros(Ng)
    
    for i in range(Ng):
        
        Vgroup = 0.2e3
        
        Vx, Vy, Vz = np.random.uniform(-1,1,3) 
        
        C = Vgroup / np.sqrt(Vx**2 + Vy**2 + Vz**2)
    
        Vel[i] = np.array([Vx * C, Vy * C, Vz * C])
        
        V[i] = LA.norm(Vel[i])
        
    return Vel
    
def PE(Pos, Mass, e, N):
    
    Pe = np.zeros(N)
    G = sc.gravitational_constant    
    
    for i in range(0,N-1):
        for j in range(i+1,N):     
            
            r = Pos[i]-Pos[j]
            m = LA.norm(r)
            
            Pe[i] += -(G*Mass[i]*Mass[j])/(m+e) # Check Pe
            Pe[j] += -(G*Mass[j]*Mass[i])/(m+e)

    return Pe

def KE(Vel, Mass, N):
    
    Ke = np.zeros(N)
    
    for i in range(0,N):
        modv = LA.norm(Vel[i])
        Ke[i] = .5*Mass[i]*modv**2

    return Ke
  
def NormV(Vel, Pos, Mass, N):
    
    Ptot = np.sum(-PE(Pos, Mass, e, N))        
    
    Ktot = np.sum(KE(Vel, Mass, N))    
    
    Tot = (2*Ktot)/Ptot
    
    l = np.random.uniform(0.95, 1)
    
    V = l*Vel/np.sqrt(Tot)
    
    return V

########################################################

def IC(Ns, Ng, R):
    
    Pos = np.zeros((Ns,3))
    V = []
    K = []
    P = []
    Mass = np.zeros(Ns)
    
    O = -1
    
    for j in range(Ng):
        i = -1
        a = int(N[j])
        
        apos = np.zeros((a,3))
        avel = np.zeros((a,3))
        amass = M(a)
        
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
                apos[i] = np.array([X,Y,Z])
                Pos[O] = apos[i] + GroupPos[j]
                
                Vx, Vy, Vz = np.random.uniform(-1,1,3)
        
                avel[i] = np.array([Vx, Vy, Vz])
                
                Mass[O] = amass[i]
                
            elif np.sqrt(X**2 + Y**2 + Z**2) > R:
               i = i - 1
               O = O - 1
               
        GroupVel = GroupV(Ng)
    
        GV = NormV(avel, apos, amass, a) + GroupVel[j]
        
        K.append(np.sum(KE(GV, amass, a)))
        P.append(np.sum(PE(apos, amass, e, a)))
        
        V.append(GV)
    
    return Pos, V, Mass, K, P
        

def IV(V, Ng, Ns):
    
    Vel = np.zeros((Ns,3))
    
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

# Dumping Data into Files

# File No.

Q = str( 1 )

np.savez('IC/IC_No'+Q+'', Pos, Vel, Mass, Ns, Ng, N, e, eta, PC, R, t_max, dt_max, t, Dump)
