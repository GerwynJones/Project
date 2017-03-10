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
Year = 365.26*(24*60*60)*(1.001)


"Defining Variables"
N = 3
t_max = Year; t = 0
dt_max = t_max/5000

v =  (2*np.pi*AU)/t_max
m = 384399000; vm = v + 1022

mass = np.array([Ms,Me,Mm])

pos = np.zeros((N,3))
vel = np.zeros((N,3))

pos[1] = np.array([0,AU,0])
pos[2] = np.array([0,AU + m,0])
vel[1] = np.array([v,0,0])
vel[2] = np.array([vm,0,0])

e = 0.0005*AU; n = 0.1

a = []; Ka = []; Pa = []
b = []; Kb = []; Pb = []
c = []; Kc = []; Pc = []

Tsum = []

T = []; dT = []

