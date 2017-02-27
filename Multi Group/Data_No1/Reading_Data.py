
from __future__ import division
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import h5py

""" Reading Data from Files """

""" Position """
with h5py.File('Position_No1.h5', 'r') as hf:
    Position = hf['Position_Data'][:]

""" Energies """
with h5py.File('Energies_No1.h5', 'r') as hf:
    Energy = hf['Energy_Data'][:]

with h5py.File('Energy_Sum_No1.h5', 'r') as hf:
    Esum = hf['Esum_Data'][:]

""" Time """
with h5py.File('Time_No1.h5', 'r') as hf:
    Time = hf['Time_Data'][:]

print Position[0]
