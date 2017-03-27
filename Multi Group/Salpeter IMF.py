# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 21:46:00 2017

@author: gezer
"""

import numpy,math
import matplotlib.pyplot as plt
import random as random
#from amuse.ic.salpeter import new_salpeter_mass_distribution

# Draw random samples from Salpeter IMF.
# N     ... number of samples.
# alpha ... power-law index.
# M_min ... lower bound of mass interval.
# M_max ... upper bound of mass interval.
def Salpeter(N, alpha_1, alpha_2, M_min, M_max):
    # Convert limits from M to logM.
    log_M_Min = math.log(M_min, 10)
    log_M_Max = math.log(M_max, 10)
    # Since Salpeter IMF decays, maximum likelihood occurs at M_min
    maxm_1 = math.pow(M_min, 1.0 - alpha_1)
    maxm_2 = math.pow(0.5, 1.0 - alpha_2)

    # Prepare array for output masses.
    MList = []

    while (len(MList) < N):
        # Draw candidate from logM interval.
        logM = random.uniform(log_M_Min,log_M_Max)
        
        M = 10**logM
        
        if M <= 0.5:
            # Compute likelihood of candidate from Salpeter SMF.
            likelihood = math.pow(M, 1.0 - alpha_1)
            # Accept randomly.
            u = random.uniform(0.0, maxm_1)
            
            if (u < likelihood):
                
                MList.append(M)
                
        if M > 0.5:
            # Compute likelihood of candidate from Salpeter SMF.
            likelihood = math.pow(M, 1.0 - alpha_2)
            # Accept randomly.
            u = random.uniform(0.0, maxm_2)
            
            if (u < likelihood):
                
                MList.append(M)
                
    Mass = numpy.array(MList)    
    
    return Mass
    
# Draw samples.
Mass = Salpeter(1e6,-5.35, 2.35, 0.1, 100.0)
# Convert to logM.
LogMass = numpy.log10(Mass)

# Plot distribution.
plt.figure(1)
plt.hist(LogMass, 200, log=True,
         range=(0.0, math.log(100.0, 10)))
         
# Overplot with Salpeter IMF.
X = []
Y = []
for n in range(101):
    logM = math.log(100.0, 10)*float(n)/100.0
    x    = 10**(logM)
    y    = 8e3*math.pow(x, 1.0 - 2.35)  # normalisation
    X.append(logM)
    Y.append(y)
plt.plot(X, Y, '-', lw=3, color='black')
plt.xlim(0.0,math.log(100.0, 10))
plt.xlabel(r'$\log M$', fontsize=20)
plt.ylabel('PDF', fontsize=20)

plt.show()