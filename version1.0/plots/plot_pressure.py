import numpy as np
import itertools
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

data = np.genfromtxt('results.csv', delimiter = ',', names = True)

rho45 = data[data['rho']==0.45]
rho65 = data[data['rho']==0.65]
rho75 = data[data['rho']==0.75]
rho85 = data[data['rho']==0.85]
rho88 = data[data['rho']==0.88]


plt.figure()
plt.xlabel(r'$\beta$', fontdict=font)
plt.ylabel(r'$\frac{\beta P}{\rho}$', fontdict=font)
plt.plot(rho45['P'], 1/rho45['T'], label=r"$\rho=0.45$")
plt.plot(rho65['P'], 1/rho65['T'], label=r"$\rho=0.65$")
plt.plot(rho75['P'], 1/rho75['T'], label=r"$\rho=0.75$")
plt.plot(rho85['P'], 1/rho85['T'], label=r"$\rho=0.85$")
plt.plot(rho88['P'], 1/rho88['T'], label=r"$\rho=0.88$")
plt.plot(data['P_verlet'], 1/data['T'], 'ro', label=r"Verlet's data")
plt.legend()
plt.show()



 


