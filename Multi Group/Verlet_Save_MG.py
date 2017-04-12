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

# File No.
Q = str(1)
# File Name
W = str('Virialised')

Position = np.zeros((Ns, len(P), 3))
Velocity = np.zeros((Ns, len(P), 3))
Acceleration = np.zeros((Ns, len(P)))
KinE = np.zeros((Ns, len(P)))
PotE = np.zeros((Ns, len(P)))
Energy = np.zeros((Ns, len(P)))
Esum = np.zeros((Ng, len(P)))
Time = np.zeros(len(P))

for j in xrange(len(P)):
    O = np.zeros(Ng+1)
    for k in xrange(Ng):
        for i in xrange(Ns):
            O[k] = O[k-1] + N[k]
            Position[i][j, :] = (P[j][i, :])/PC
            Velocity[i][j, :] = V[j][i, :]
            Acceleration[i, j] = A[j][i]
            KinE[i, j] = KE[j][i]
            PotE[i, j] = PE[j][i]
            Energy[i, j] = TE[j][i]
            Esum[k, j] = np.sum(TE[j][O[k-1]:O[k]])
            Time[j] = T[j]/Year

# Dumping Data into Files

# Position
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Position_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Position_Data",  data=Position)

# Velocity
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Velocity_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Velocity_Data",  data=Velocity)

# Acceleration
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Acceleration_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Acceleration_Data",  data=Acceleration)

# Kinetic Energy
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Kin_Energy_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Kin_Energy_Data",  data=KinE)

# Potential Energy
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Pot_Energy_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Pot_Energy_Data",  data=PotE)

# Energy
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Energy_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Energy_Data",  data=Energy)

# Total Energy
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Energy_Sum_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Esum_Data",  data=Esum)

# Time
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Time_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Time_Data",  data=Time)

# Mass
with h5py.File('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Mass_No'+Q+'.h5', 'w') as hf:
    hf.create_dataset("Mass_Data",  data=Mass)


# Percentage of Virial

Percent = (Type/Virial)*100

print(Percent)

String = np.array(['Percent = '])

Load = np.column_stack((String, Percent))

np.savetxt('/media/gerwyn/Linux Files/Work/PLF/'+W+'/Data_No'+Q+'/Percentage.txt', Load, delimiter=" ", fmt="%s")
