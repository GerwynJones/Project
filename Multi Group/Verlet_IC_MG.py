#!/usr/bin/env python
"""
Created on Thu Nov 24 22:29:44 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
import math
import scipy.constants as sc
from numpy import linalg as LA
# import h5py

###############################################

# Defining Variables
# Size
AU = sc.astronomical_unit
PC = 206265*AU
R = 200*AU

# No.Of.groups
Ng = int(1e3)  #10

# Dumping Number
Dump = 25

# Duration
Year = sc.Julian_year
t_max = 1e6*Year; t = 0; dt_max = 1.5*Year

# Initial Conditions
# Constants
e = 0.05*AU; eta = 5

############################################

def Salpeter(N, alpha, M_min, M_max):
    # Convert limits from M to logM.
    log_M_Min = math.log(M_min, 10)
    log_M_Max = math.log(M_max, 10)

    C = 1.35

    # Since Salpeter IMF decays, maximum likelihood occurs at M_min
    maxm = math.pow(M_min, 1.0 - alpha)/C

    # Prepare array for output masses.
    MList = []

    while (len(MList) < N):
        # Draw candidate from logM interval.
        logM = np.random.uniform(log_M_Min,log_M_Max)

        M = 10**logM

        # Compute likelihood of candidate from Salpeter SMF.
        likelihood = math.pow(M, 1.0 - alpha)
        # Random
        u = np.random.uniform(0.0, maxm)

        if (u < likelihood):

            MList.append(M)

    Mass = np.array(MList)

    return Mass

def Kroupa(N, alpha, M_min, M_max):
    # Convert limits from M to logM.
    log_M_Min = math.log(M_min, 10)
    log_M_Max = math.log(M_max, 10)

    alpha_1 = alpha[0]
    alpha_2 = alpha[1]

    C_1 = 0.35
    C_2 = 1.35

    # Since Kroupa IMF decays, maximum likelihood occurs at M_min
    maxm_1 = math.pow(M_min, 1.0 - alpha_1)/C_1
    maxm_2 = math.pow(M_min + 0.02, 1.0 - alpha_2)/C_2

    # Prepare array for output masses.
    MList = []

    while (len(MList) < N):
        # Draw candidate from logM interval.
        logM = np.random.uniform(log_M_Min,log_M_Max)

        M = 10**logM

        if M <= 0.5:
            # Compute likelihood of candidate from Kroupa IMF.
            likelihood = math.pow(M, 1.0 - alpha_1)
            # Random
            u = np.random.uniform(0.0, maxm_1)

            if (u < likelihood):

                MList.append(M)

        if M > 0.5:
            # Compute likelihood of candidate from Kroupa IMF.
            likelihood = math.pow(M, 1.0 - alpha_2)
            # Random
            u = np.random.uniform(0.0, maxm_2)

            if (u < likelihood):

                MList.append(M)

    Mass = np.array(MList)

    return Mass

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
            
            S = 200  #np.random.randint(3, 6)
            
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

            S = 200  #np.random.randint(3, 6)
            
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
    G = sc.gravitational_constant    
    
    for i in xrange(0, N - 1):
        for j in xrange(i + 1, N):     
            
            r = Pos[i]-Pos[j]
            m = LA.norm(r)
            
            Pe[i] += -(G*Mass[i]*Mass[j])/(m+e) # Check Pe
#            Pe[j] += -(G*Mass[j]*Mass[i])/(m+e)

    return Pe

def KE(Vel, Mass, N):
    
    Ke = np.zeros(N)
    
    for i in xrange(0, N):
        modv = LA.norm(Vel[i])
        Ke[i] = .5*Mass[i]*modv**2

    return Ke
  
def NormV(Vel, Pos, Mass, N):
    
    Ptot = np.sum(-PE(Pos, Mass, e, N))        
    
    " This is not Virial as it doesnt take into account the applied group velocity "
    
    Ktot = np.sum(KE(Vel, Mass, N))
    
    Tot = (2*Ktot)/Ptot
    
    V = Vel/np.sqrt(Tot)
    
    return V

########################################################

def IC(Ns, Ng, R, GroupPos, N, Method):
    """ Creating the initial conditions with random points within a cylinder """    
    
    Pos = np.zeros((Ns, 3))
    V = []
    K = []
    P = []
    Mass = np.zeros(Ns)
    
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
    
        GV = NormV(avel, apos, amass, a) + GroupVel[j]
        
        Velf = NormV(GV, apos, amass, a)
        
        K.append(np.sum(KE(Velf, amass, a)))
        P.append(np.sum(PE(apos, amass, e, a)))
        
        V.append(Velf)
    
    return Pos, V, Mass, K, P
        

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

Pos, V, Mass, KinE, PotE = IC(Ns, Ng, R, GroupPos, N, Kroupa)

Vel = IV(V, Ng, Ns, N)



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(Pos[:, 0]/PC, Pos[:, 1]/PC, Pos[:, 2]/PC, color='green', linestyle='None', marker='.')

ax.set_xlabel("Distance (Pc)", fontsize=14)
ax.set_ylabel("Distance (Pc)", fontsize=14)
ax.set_zlabel("Distance (Pc)", fontsize=14)
plt.legend(loc='best')

fig = plt.figure()

plt.plot(Pos[:, 0]/PC, Pos[:, 1]/PC, color='green', linestyle='None', marker='.', markersize=9)

plt.xlabel("Distance (Pc)", fontsize=18)
plt.ylabel("Distance (Pc)", fontsize=18)
plt.legend(loc='best')
plt.savefig('Graphs/Graph of IC Cylinder_BOLD.png', bbox_inches='tight')


fig = plt.figure()

plt.plot(Pos[:, 1]/PC, Pos[:, 2]/PC, color='green', linestyle='None', marker='.', markersize=9)

plt.xlabel("Distance (Pc)", fontsize=18)
plt.ylabel("Distance (Pc)", fontsize=18)
plt.legend(loc='best')
plt.savefig('Graphs/Graph of IC Circle_BOLD.png', bbox_inches='tight')

plt.show()




# Dumping Data into Files

# File No.

# Q = str( 1 )
#
# # Position
# with h5py.File('IC_No'+Q+'/Position_No'+Q+'.h5', 'w') as hf:
#    hf.create_dataset("Position_Data",  data=Pos)
#
# # Velocity
# with h5py.File('IC_No'+Q+'/Velocity_No'+Q+'.h5', 'w') as hf:
#    hf.create_dataset("Velocity_Data",  data=Vel)
#    
# # Mass
# with h5py.File('IC_No'+Q+'/Mass_No'+Q+'.h5', 'w') as hf:
#   hf.create_dataset("Mass_Data",  data=Mass)
#
# # N
# with h5py.File('IC_No'+Q+'/N_No'+Q+'.h5', 'w') as hf:
#   hf.create_dataset("N_Data",  data=N)
#
# # NGroup
# with open('IC_No'+Q+'/Ng_No'+Q+'.txt', 'w') as f:
#  f.write('%d' % Ng)
#
# # Ns
# with open('IC_No'+Q+'/Ns_No'+Q+'.txt', 'w') as f:
#  f.write('%d' % Ns)
#
# # Dump
# with open('IC_No'+Q+'/Dump_No'+Q+'.txt', 'w') as f:
#  f.write('%d' % Dump)
#
# # T_max
# with open('IC_No'+Q+'/Tmax_No'+Q+'.txt', 'w') as f:
#  f.write('%d' % t_max)
#
# # T
# with open('IC_No'+Q+'/t_No'+Q+'.txt', 'w') as f:
#  f.write('%d' % t)
#
# # dT
# with open('IC_No'+Q+'/dT_No'+Q+'.txt', 'w') as f:
#  f.write('%d' % dt_max)
#
# # eta
# with open('IC_No'+Q+'/eta_No'+Q+'.txt', 'w') as f:
#  f.write('%f' % eta)
#
# # e
# with open('IC_No'+Q+'/e_No'+Q+'.txt', 'w') as f:
#  f.write('%d' % e)



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
