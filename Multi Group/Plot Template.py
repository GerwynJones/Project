#!/usr/bin/env python
"""
Created on Thu Mar  2 23:50:10 2017

@author: Admin
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plt.plot(Pos[:, 0], Pos[:, 1], Pos[:, 2], color='green', linestyle='None', marker='.')

plt.xlabel("")
plt.legend(loc='best')


fig = plt.figure()

plt.plot(Pos[:, 0]/PC, Pos[:, 1]/PC, color='green', linestyle='None', marker='.')

plt.xlabel("Distance (Pc)")
plt.ylabel("Distance (Pc)")
plt.legend(loc='best')
plt.savefig('Graphs/Graph of IC Cylinder_2.png', bbox_inches='tight')


fig = plt.figure()

plt.plot(Pos[:, 1]/PC, Pos[:, 2]/PC, color='green', linestyle='None', marker='.')

plt.xlabel("Distance (Pc)")
plt.ylabel("Distance (Pc)")
plt.legend(loc='best')
plt.savefig('Graphs/Graph of IC Circle_2.png', bbox_inches='tight')

plt.show()
