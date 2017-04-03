# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 22:29:44 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
import scipy.constants as sc

AU = sc.astronomical_unit
Ms = 1.989e30
Me = 5.972e24
Mm = 7.342e22
Year = sc.Julian_year


"Defining Variables"
N = 3
t_max = 10*Year; t = 0
dt_max = Year/50

v = (2*np.pi*AU)/Year; vm = v + 1022
m = 384399000

mass = np.array([Ms, Me, Mm])

pos = np.zeros((N, 3))
vel = np.zeros((N, 3))

pos[1] = np.array([0, AU, 0])
vel[1] = np.array([v, 0, 0])
pos[2] = np.array([0, AU + m, 0])
vel[2] = np.array([vm, 0, 0])

e = 0.0005*AU; n = 0.1

a = []; b = []

ea = []; eb = []

c = []; ec = []

Tsum = []; T = []

