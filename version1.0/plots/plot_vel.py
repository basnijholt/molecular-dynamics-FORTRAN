import numpy as np
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
#        'size'   : 16,
        }

vel_last = np.loadtxt('vel.dat') 
vel_IC = np.loadtxt('vel_IC.dat')

v1 = np.sqrt((vel_IC**2).sum(axis=1))
v2 = np.sqrt((vel_last**2).sum(axis=1))

bins = 40

plt.figure()
plt.subplot(2, 1, 1)
plt.ylabel(r'Frequency', fontdict=font)
#plt.xlabel(r'Velocity', fontdict=font)

count, bins, ignored = plt.hist(v1, bins, normed=True, label=r'Initial velocities')
plt.plot(bins, 1/(np.std(v1) * np.sqrt(2 * np.pi))* np.exp( - (bins - np.mean(v1))**2 / (2 *np.std(v1)**2) ), linewidth=2, color='r')
plt.legend()

plt.subplot(2, 1, 2)
plt.ylabel(r'Frequency', fontdict=font)
plt.xlabel(r'Velocity ($\sigma / \Delta t$)', fontdict=font)
count, bins, ignored = plt.hist(v2, bins, normed=True, label=r'Velocities at $t=2000$')
#plt.plot(bins, 1/(np.std(v2) * np.sqrt(2 * np.pi))* np.exp( - (bins - np.mean(v2))**2 / (2 *np.std(v2)**2) ), linewidth=2, color='r')
plt.legend()

plt.show()
