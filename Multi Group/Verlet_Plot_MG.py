#!/usr/bin/env python
"""
Created on Thu Nov 24 22:59:34 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
import scipy.constants as sc
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import h5py

from Verlet_main_MG import *

plt.close('all')

print("Starting Plot File")

Position = np.zeros((Ns,len(P),3))
Energy = np.zeros((Ns,len(P)))
Esum = np.zeros((Ng,len(P)))
Time = np.zeros(len(P))

T_max = np.max(dT) 
T_min = np.min(dT)
T_Ratio = T_max/T_min

for j in range(len(P)):
    O = np.zeros(Ng+1)
    for k in range(Ng): 
        for i in range(Ns):
            O[k] = O[k-1] + N[k]
            Position[i][j,:] = P[j][i,:]
            Energy[i,j] = E[j][i]
            Esum[k,j] = np.sum(E[j][O[k-1]:O[k]])
            Time[j] = T[j]
    
        
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')          
        
for i in range(Ns):    
    plt.plot(Position[i][:,0],Position[i][:,1],Position[i][:,2]) 

plt.xlabel("x")
plt.legend(loc = 'best')

plt.figure()

for i in range(Ns):    
    j = i + 1
    plt.plot(Position[i][:,0],Position[i][:,1])
    
plt.ylabel(r'Distance $(AU)$')
plt.xlabel(r'Distance $(AU)$')
#plt.savefig('Graphs/Graph of 2D MG.png', bbox_inches='tight')
plt.legend(loc = 'best')


O = np.zeros(Ng + 1)

for k in range(Ng):
    plt.figure()
    for i in range(int(N[k])):
        O[k] = O[k-1] + N[k]
        i = i + O[k-1] 
        
        j = i + 1
        
        plt.plot(Time, Energy[i,:]/Esum[k,:], label = 'star %s' % j)
        #  /Esum[k,:]
    plt.xlabel("Time (yrs)")
    plt.legend(loc = 'best')
    
#plt.savefig('Graphs/Graph of Energy MG.png', bbox_inches='tight')


plt.plot(Time, E)
plt.xlabel("Time (yrs)")
plt.ylabel("Energy (J)")
plt.title("Graph of Total energy of system")
plt.show()

# Dumping Data into Files

# File No.
Q = str(2)

# Position
with h5py.File('Data_No'+Q+'/Position_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Position_Data",  data=Position)

# Energies
with h5py.File('Data_No'+Q+'/Energies_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Energy_Data",  data=Energy)

with h5py.File('Data_No'+Q+'/Energy_Sum_No'+Q+'.h5', 'w') as hf:
   hf.create_dataset("Esum_Data",  data=Esum)

# Time
with h5py.File('Data_No'+Q+'/Time_No'+Q+'.h5', 'w') as hf:
   hf.create_dataset("Time_Data",  data=Time)

# Mass
with h5py.File('Data_No'+Q+'/Mass_No'+Q+'.h5', 'w') as hf:
   hf.create_dataset("Mass_Data",  data=Mass)


Year = sc.Julian_year

TR = np.array([(T_min/Year),(T_max/Year),T_Ratio])
ST = np.array(['T Min','T Max','T Ratio'])

LoadTime = np.column_stack((ST, TR))

np.savetxt('Data_No'+Q+'/Time.txt', LoadTime, delimiter=" ", fmt="%s")

