# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 22:29:44 2016

@author: Admin
"""
from __future__ import division
import numpy as np 
import scipy.constants as sc

AU = sc.astronomical_unit; A = 200*AU
Ms = 1.989e30
Me = 5.972e24
Year = 365.26*(24*60*60)


"Defining Variables"
N = 2
t_max = Year*(1.25e3)  #1.25e3 for dt_grav and .65e3 for dt_max
dt_max = t_max/200

V =  (2*np.pi*A)/(Year*(5.5e3))

mass = np.array([Ms,Me])

pos = np.zeros((N,3))
vel = np.zeros((N,3))

pos[1] = np.array([0,A,0])
vel[1] = np.array([V,0,0])

e = 0.05*AU; n = 0.05

a = []; Ta = []
b = []; Tb = []

Tsum = []

T = []; dT = []