#!/usr/bin/env python
"""
Created on Fri Mar 24 21:46:00 2017

@author: Admin
"""
from __future__ import division
import numpy, math
import matplotlib.pyplot as plt
import random as random

# Draw random samples from Kroupa IMF.
# N     = number of samples.
# alpha = power-law index.
# M_min = lower bound of mass interval.
# M_max = upper bound of mass interval.
def Kroupa(N, alpha, M_min, M_max):
    # Convert limits from M to logM.
    log_M_Min = math.log(M_min, 10)
    log_M_Max = math.log(M_max, 10)

    alpha_1 = alpha[0]
    alpha_2 = alpha[1]

    C_1 = 0.3
    C_2 = 1.3

    # Since Kroupa IMF decays, maximum likelihood occurs at M_min
    maxm_1 = math.pow(M_min, 1.0 - alpha_1)/C_1

    maxm_2 = math.pow(M_min + 0.01, 1.0 - alpha_2)/C_2

    # Prepare array for output masses.
    MList = []

    while len(MList) < N:
        # Draw candidate from logM interval.
        logM = random.uniform(log_M_Min,log_M_Max)

        M = 10**logM

        if M <= 0.5:
            # Compute likelihood of candidate from Kroupa IMF.
            likelihood = math.pow(M, 1.0 - alpha_1)
            # Random
            u = random.uniform(0.0, maxm_1)

            if (u < likelihood):

                MList.append(M)

        if M > 0.5:
            # Compute likelihood of candidate from Kroupa IMF.
            likelihood = math.pow(M, 1.0 - alpha_2)
            # Random
            u = random.uniform(0.0, maxm_2)

            if u < likelihood:

                MList.append(M)

    Mass = numpy.array(MList)

    return Mass

"""
# Draw samples from Kroupa
Mass = Kroupa(1e6, numpy.array([1.3, 2.3]), 0.1, 100.0)
# Convert to logM.
LogMass = numpy.log10(Mass)

# Plot distribution.
plt.figure(1)
plt.hist(LogMass, 100, log=True,
         range=(math.log(0.1, 10), math.log(100.0, 10)))

# Overplot with Kroupa IMF.
X = []
Y = []

logM = numpy.linspace(math.log(0.1, 10), math.log(100, 10), 100)

for n in range(len(logM)):
    C_1 = 0.3
    C_2 = 1.3

    if logM[n] <= math.log(0.5, 10):

        x    = 9.5**(logM[n])
        y    = 7e3*math.pow(x, 1.0 - 1.3)/C_1  # normalisation
        X.append(logM[n])
        Y.append(y)

    if logM[n] > math.log(0.5, 10):

        x    = 10**(logM[n])
        y    = 15.5e3*math.pow(x, 1.0 - 2.3)/C_2  # normalisation
        X.append(logM[n])
        Y.append(y)

plt.plot(X, Y, '-', lw=5, color='red')
plt.xlim(math.log(0.1, 10),math.log(100.0, 10))
plt.xlabel(r'$\log M$', fontsize=18)
plt.ylabel('PDF', fontsize=18)

plt.savefig('Initial conditions/Kroupa_IMF.png', bbox_inches='tight')

plt.show()
"""
