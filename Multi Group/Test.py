#!/usr/bin/env python
"""
Created on Thu Mar  2 23:50:10 2017

@author: Admin
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from mpl_toolkits.mplot3d import Axes3D

plt.close('all')

# Dumping Data into Files
# File No.
Q = str(6)

Type = 0.4
Virial = 1/2

# Percentage of Virial

Percent = (Type/Virial)*100

print(Percent)

String = np.array(['Percent = '])

Load = np.column_stack((String, Percent))

np.savetxt('/media/gerwyn/Linux Files/Work/PLF/Data_No'+Q+'/Percentage.txt', Load, delimiter=" ", fmt="%s")
