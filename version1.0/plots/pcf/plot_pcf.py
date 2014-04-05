import numpy as np
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

solid = np.loadtxt('T0.01rho0.85.dat') 
liquid = np.loadtxt('T0.30rho0.85.dat')
gas = np.loadtxt('T1.0rho0.85.dat')

r_gas = gas[:,0]
g_gas = gas[:,1]
r_liquid = liquid[:,0]
g_liquid = liquid[:,1]
r_solid = solid[:,0]
g_solid = solid[:,1]

plt.figure()

plt.subplot(3, 1, 1)
plt.title(r'Pair correlation function $g(r)$', fontdict=font)
plt.ylabel(r'$g(r)$', fontdict=font)
#plt.xlabel(r'$r/\sigma$', fontdict=font)
plt.hlines(1, 0, 4.5, linestyles='dashed')
plt.plot(r_gas, g_gas, 'r', label=r'$g(r)$ for a gas')
plt.legend()

plt.subplot(3, 1, 2)
plt.ylabel(r'$g(r)$', fontdict=font)
#plt.xlabel(r'$r/\sigma$', fontdict=font)
plt.hlines(1, 0, 4.5, linestyles='dashed')
plt.plot(r_liquid, g_liquid, 'b', label=r'$g(r)$ for a liquid')
plt.legend()

plt.subplot(3, 1, 3)
plt.ylabel(r'$g(r)$', fontdict=font)
plt.xlabel(r'$r/\sigma$', fontdict=font)
plt.hlines(1, 0, 4.5, linestyles='dashed')
plt.plot(r_solid, g_solid, 'y', label=r'$g(r)$ for a solid')
plt.legend()
plt.show()
