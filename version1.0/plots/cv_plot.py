import numpy as np
import itertools
import matplotlib.pyplot as plt

data = np.genfromtxt('results.csv', delimiter = ',', names = True)

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

# data
rho45 = data[data['rho']==0.45]
rho65 = data[data['rho']==0.65]
rho75 = data[data['rho']==0.75]
rho85 = data[(data['rho']==0.85)&(data['T']<=2.0)]
rho88 = data[data['rho']==0.88]



plt.figure()
plt.ylabel(r'Specific heat $C_V$', fontdict=font)
plt.xlabel(r'Temperature ($\epsilon / k$)', fontdict=font)

# extra punten berekend, de data niet verder opgeslagen
plt.plot([7.981, 6, 8.5], [1319.88, 1323.15, 1317.73], 'mo', label=r"$\rho=0.1$")

plt.plot(rho45['T'], rho45['cv'], 'bo', label=r"$\rho=0.45$")
plt.plot(rho85['T'], rho85['cv'], 'yo', label=r"$\rho=0.85$")
plt.plot(rho88['T'], rho88['cv'], 'ro', label=r"$\rho=0.88$")

# nog wat extra punten met lage T
plt.plot([0.25, 0.05, 0.1, 0.06, 0.05, 0.055], [2304, 2504, 2399, 2429, 2614, 2680], 'ro')

#lijnen om de ideale Cv's aan te geven
plt.axhline(y=3*864, xmin=0.01, xmax=0.3, color='y', marker='D',linewidth=2, label=r"$\frac{3}{2}Nk$")
plt.axhline(y=3*864/2, xmin=0.7, xmax=0.95, color='y', marker='o',linewidth=2, label=r"$3Nk$")


plt.legend()
plt.show()



 


