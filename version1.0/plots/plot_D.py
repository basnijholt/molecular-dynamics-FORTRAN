import numpy as np
import itertools
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
#        'size'   : 16,
        }

data = np.genfromtxt('results.csv', delimiter = ',', names = True)

rho45 = data[data['rho']==0.45]
rho65 = data[data['rho']==0.65]
rho75 = data[data['rho']==0.75]
rho85 = data[data['rho']==0.85]
rho88 = data[data['rho']==0.88]


plt.figure()
plt.subplot(2,1,1)
plt.xlabel(r'Pressure ($\sigma^3 / \epsilon$)', fontdict=font)
plt.ylabel(r'Diffusion coefficient ($\sigma^2/\Delta t$)', fontdict=font)
plt.plot(rho85['P'], 	rho85['D'], 'ro-')
plt.subplot(2,1,2)
plt.xlabel(r'Temperature ($\epsilon / k$)', fontdict=font)
plt.ylabel(r'D coefficient ($\sigma^2/\Delta t$)', fontdict=font)
plt.plot(rho85['T'], 	rho85['D'], 'ro-')
plt.legend()
plt.show()



 


