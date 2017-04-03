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

# Dumping Data into Files
# File No.
Q = str(6)
plt.close('all')

print("Starting Plot File")

# Size
AU = sc.astronomical_unit
PC = 206265*AU

Position = np.zeros((Ns, len(P), 3))
Energy = np.zeros((Ns, len(P)))
Esum = np.zeros((Ng, len(P)))
Time = np.zeros(len(P))

for j in xrange(len(P)):
    O = np.zeros(Ng+1)
    for k in xrange(Ng):
        for i in xrange(Ns):
            O[k] = O[k-1] + N[k]
            Position[i][j, :] = (P[j][i, :])/PC
            Energy[i, j] = E[j][i]
            Esum[k, j] = np.sum(E[j][O[k-1]:O[k]])
            Time[j] = T[j]/Year
    
        
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')          
        
for i in xrange(Ns):
    plt.plot(Position[i][:, 0], Position[i][:, 1], Position[i][:, 2])

plt.xlabel("x")
plt.legend(loc='best')


plt.figure()

for i in xrange(Ns):
    j = i + 1
    plt.plot(Position[i][:, 0], Position[i][:, 1])
    
plt.ylabel(r'Distance $(Pc)$')
plt.xlabel(r'Distance $(Pc)$')
plt.savefig('Graphs/Graph of 2D MG NO_'+Q+'.png', bbox_inches='tight')
plt.legend(loc='best')


O = np.zeros(Ng + 1)
for k in xrange(Ng):
    plt.figure()
    l = str(k + 1)
    for i in xrange(int(N[k])):
        O[k] = O[k-1] + N[k]
        i = i + O[k-1]
        j = i + 1
        
        plt.plot(Time, Energy[i, :], label='star %s' % j)
    plt.ylabel("Energy (J)")
    plt.xlabel("Time (yrs)")
    plt.legend(loc='best')
    plt.savefig('Graphs/Graph of Energy MG for group '+l+' NO_'+Q+'.png', bbox_inches='tight')

    plt.figure()
    plt.plot(Time, Esum[k, :])
    plt.xlabel("Time (yrs)")
    plt.ylabel("Energy (J)")
    plt.savefig('Graphs/Graph of Total energy of system for group '+l+' NO_'+Q+'.png', bbox_inches='tight')

plt.show()

# Position
with h5py.File('/home/gerwyn/Documents/Project-Large-Files/Data_No'+Q+'/Position_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Position_Data",  data=Position)

# Energies
with h5py.File('/home/gerwyn/Documents/Project-Large-Files/Data_No'+Q+'/Energies_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Energy_Data",  data=Energy)

with h5py.File('/home/gerwyn/Documents/Project-Large-Files/Data_No'+Q+'/Energy_Sum_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Esum_Data",  data=Esum)

# Time
with h5py.File('/home/gerwyn/Documents/Project-Large-Files/Data_No'+Q+'/Time_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Time_Data",  data=Time)

# Mass
with h5py.File('/home/gerwyn/Documents/Project-Large-Files/Data_No'+Q+'/Mass_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Mass_Data",  data=Mass)

Year = sc.Julian_year

T_max = np.max(dT)
T_min = np.min(dT)
T_Ratio = T_max/T_min

TR = np.array([(T_min/Year), (T_max/Year), T_Ratio])
ST = np.array(['T Min', 'T Max', 'T Ratio'])

LoadTime = np.column_stack((ST, TR))

np.savetxt('/home/gerwyn/Documents/Project-Large-Files/Data_No'+Q+'/Time.txt', LoadTime, delimiter=" ", fmt="%s")

