import numpy as np
import matplotlib.pyplot as plt
font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

T = np.loadtxt('t.dat') 


U = np.loadtxt('u.dat')

time = np.linspace(0,6000,6000)

plt.figure()

plt.subplot(2, 1, 1)
#for col_name in data.dtype.names:
plt.gca().axes.get_xaxis().set_visible(False)
plt.ylabel(r'Total energy $(\epsilon)$', fontdict=font)

plt.plot(time[0:1000], U[0:1000], 'r')


# stop temperature rescaling
plt.annotate(
    '', xy=(600, -4200), xycoords = 'data',
    xytext = (500, -3871), textcoords = 'data',
    arrowprops = {'arrowstyle':'<-'})

plt.annotate(
    r'$t_{eq}$', xy=(605, -4270), xycoords = 'data',
    xytext = (0, 0), textcoords = 'offset points', size=16)



plt.subplot(2, 1, 2)
plt.ylim(0,1.1)
plt.xlabel(r'time ($\Delta t$)', fontdict=font)
plt.ylabel(r'Temperature $(\epsilon / k)$', fontdict=font)

plt.plot(time[0:1000], T[0:1000], 'g')


# stop temperature rescaling
plt.annotate(
    '', xy=(600, 0.7), xycoords = 'data',
    xytext = (500, 0.98), textcoords = 'data',
    arrowprops = {'arrowstyle':'<-'})

plt.annotate(
    r'$t_{eq}$', xy=(605, 0.65), xycoords = 'data',
    xytext = (0, 0), textcoords = 'offset points', size=16)


plt.legend()
plt.show()



 
