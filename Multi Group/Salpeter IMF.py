# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 21:46:00 2017

@author: Admin
"""
from __future__ import division
import numpy,math
import matplotlib.pyplot as plt
import random as random

# Draw random samples from Salpeter IMF.
# N     ... number of samples.
# alpha ... power-law index.
# M_min ... lower bound of mass interval.
# M_max ... upper bound of mass interval.
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
        logM = random.uniform(log_M_Min,log_M_Max)
        
        M = 10**logM

        # Compute likelihood of candidate from Salpeter SMF.
        likelihood = math.pow(M, 1.0 - alpha)
        # Random
        u = random.uniform(0.0, maxm)

        if (u < likelihood):

            MList.append(M)

    Mass = numpy.array(MList)    
    
    return Mass
    
# Draw samples.
Mass = Salpeter(1e6, 2.35, 0.4, 100.0)
# Convert to logM.
LogMass = numpy.log10(Mass)

# Plot distribution.
plt.figure(1)
plt.hist(LogMass, 200, log=True,
         range=(math.log(0.4, 10), math.log(100.0, 10)))
         
# Overplot with Salpeter IMF.

X = []
Y = []

logM = numpy.linspace(math.log(0.4, 10), math.log(100, 10),100)

for n in range(len(logM)):
    C = 1.35

    x    = 10**(logM[n])
    y    = 1.5e4*math.pow(x, 1.0 - 2.35)/C  # normalisation
    X.append(logM[n])
    Y.append(y)

plt.plot(X, Y, '-', lw=3, color='red')
plt.xlim(math.log(0.5, 10), math.log(100.0, 10))
plt.xlabel(r'$\log M$', fontsize=15)
plt.ylabel('PDF', fontsize=15)

plt.savefig('Graphs/Salpeter_IMF.png', bbox_inches='tight')

plt.show()
