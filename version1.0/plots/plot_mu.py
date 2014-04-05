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
plt.xlabel(r'Temperature ($\epsilon / k$)', fontdict=font)
plt.ylabel(r'Diffusion coefficient ($\sigma^2/\Delta t$)', fontdict=font)
#plt.xlim(-1.6,0.5)
#plt.plot(rho45['mu'], 1/rho45['T'], label=r"$\rho=0.45$")
#plt.plot(rho65['mu'], 1/rho65['T'], label=r"$\rho=0.65$")
#plt.plot(rho75['mu'], 1/rho75['T'], label=r"$\rho=0.75$")
plt.plot(rho85['T'][2:], 	rho85['mu'][2:], 'ro-')
plt.xscale('symlog')
plt.yscale('symlog')
plt.ylim(-.0005,0)
#plt.errorbar(rho85['mu'], 1/rho85['T'], yerr=10*rho85['mu_error'], linestyle="None", marker="None", color="green")
#plt.plot(rho88['mu'], 1/rho88['T'], label=r"$\rho=0.88$")
plt.legend()
plt.show()



 


