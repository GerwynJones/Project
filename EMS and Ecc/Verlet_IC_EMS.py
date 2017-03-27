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
Year = sc.Julian_year



"Defining Variables"
N = 2
t_max = 100*Year; t = 0
dt_max = Year/5

v = 29754.7

mass = np.array([Ms,Me])

pos = np.zeros((N,3))
vel = np.zeros((N,3))

pos[1] = np.array([0,AU,0])
vel[1] = np.array([v,0,0])

e = 0.005*AU; n = 0.1

a = []; b = []

ea = []; eb = []

Tsum = []; T = []

N1 = []; N2 = []
