import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


fcc = np.loadtxt('fcc.dat')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n = 100
xs = fcc[:,0]
ys = fcc[:,1]
zs = fcc[:,2]
ax.scatter(xs, ys, zs, c='r', marker='o')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

