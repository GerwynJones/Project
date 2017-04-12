#!/usr/bin/env python
"""
Created on Thu Mar  2 23:50:10 2017

@author: Admin
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


# Plotting Graphs
plt.close('all')

print("Starting Plot File")

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

for i in xrange(Ns):
    plt.plot(Position[i][:, 0], Position[i][:, 1], Position[i][:, 2])

ax.set_xlabel("Distance (Pc)", fontsize=14)
ax.set_ylabel("Distance (Pc)", fontsize=14)
ax.set_zlabel("Distance (Pc)", fontsize=14)
plt.legend(loc='best')

plt.figure()

for i in xrange(Ns):
    j = i + 1
    plt.plot(Position[i][:, 0], Position[i][:, 1])

plt.ylabel(r'Distance $(Pc)$', fontsize=14)
plt.xlabel(r'Distance $(Pc)$', fontsize=14)
plt.savefig('Graphs/Graph of 2D MG NO_'+Q+'.png', bbox_inches='tight')
plt.legend(loc='best')


O = np.zeros(Ng + 1)
for k in xrange(Ng):
    plt.figure()
    l = str(k + 1)
    for i in xrange(int(N[k])):
        O[k] = O[k-1] + N[k]
        i = i + O[k-1]
        j = i + 1

        plt.plot(Time, Energy[i, :], label='star %s' % j)

    plt.ylabel("Energy (J)", fontsize=14)
    plt.xlabel("Time (yrs)", fontsize=14)
    plt.legend(loc='best')
    plt.savefig('Graphs/Graph of Energy MG for group '+l+' NO_'+Q+'.png', bbox_inches='tight')

    plt.figure()

    plt.plot(Time, Esum[k, :])

    plt.xlabel("Time (yrs)", fontsize=14)
    plt.ylabel("Energy (J)", fontsize=14)
    plt.savefig('Graphs/Graph of Total energy of system for group '+l+' NO_'+Q+'.png', bbox_inches='tight')

plt.show()

