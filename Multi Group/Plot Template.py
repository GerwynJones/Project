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

ax.plot(Pos[:, 0]/PC, Pos[:, 1]/PC, Pos[:, 2]/PC, color='green', linestyle='None', marker='.')

ax.set_xlabel("Distance (Pc)", fontsize=14)
ax.set_ylabel("Distance (Pc)", fontsize=14)
ax.set_zlabel("Distance (Pc)", fontsize=14)
ax.set_xlim(-0.05, 1.05)
plt.legend(loc='best')


fig = plt.figure()

plt.plot(Pos[:, 0]/PC, Pos[:, 1]/PC, color='green', linestyle='None', marker='.', markersize=9)

plt.xlabel("Distance (Pc)", fontsize=18)
plt.ylabel("Distance (Pc)", fontsize=18)
plt.xlim(-0.05, 1.05)
plt.legend(loc='best')
plt.savefig('Graphs/Graph of IC Cylinder_BOLD.png', bbox_inches='tight')


fig = plt.figure()

plt.plot(Pos[:, 1]/PC, Pos[:, 2]/PC, color='green', linestyle='None', marker='.', markersize=9)

plt.xlabel("Distance (Pc)", fontsize=18)
plt.ylabel("Distance (Pc)", fontsize=18)
plt.legend(loc='best')
plt.savefig('Graphs/Graph of IC Circle_BOLD.png', bbox_inches='tight')

plt.show()
