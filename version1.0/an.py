import numpy as np
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

array = np.loadtxt('array.dat') #T, kin_energy, U, PkT, disp_sq
T = array[:,0]
K = array[:,1]
U = array[:,2]
P = array[:,3]
disp = array[:,4]

data = np.loadtxt('data.dat') #rho, T_IC, eq_time, time_end, dt, r_cut, r_verlet, N_bin
rho = data[0]
T_IC = data[1]
eq_time = data[2]
time_end = data[3]
dt = data[4]
r_cut = data[5]
r_verlet = data[6]
N_bin = data[7]
N = data[8]

entries = np.size(T)

# datablock!
data_block_size = 50
def data_block(x):
 correct_length = entries-entries%data_block_size # make sure the list is an integer multiple of data_block_size
 xx = x[:correct_length].reshape(-1,data_block_size)
 avgs = np.mean(xx, axis = 0)
 var = np.std(xx,axis = 0)
 avg = np.mean(avgs, axis = 0)
 error = np.std(avgs)/np.sqrt(np.size(avgs))
 return avgs, var, avg, error

T_avgs, T_var, T_avg, T_error = data_block(T)
K_avgs, K_var, K_avg, K_error = data_block(K)
U_avgs, U_var, U_avg, U_error = data_block(U)
P_avgs, P_var, P_avg, P_error = data_block(P)

cv = 1.0/(2.0/(3.0*864.0)-(K_var/K_avgs)**2)
cv_avg = np.mean(cv)
cv_error = np.std(cv)/np.sqrt(np.size(cv))

# load mu
mu = np.loadtxt('mu.dat')
mu_avg = np.mean(mu, axis = 0)
mu_error = np.std(mu)/np.sqrt(np.size(mu))

# de least-square methode voor de fit van de diffusie coefficient
n = np.size(disp)
x = np.arange(n)
A = np.vstack([x, np.ones(n)]).T
model, resid = np.linalg.lstsq(A, disp)[:2]
r2 = 1 - resid / (np.size(disp) * np.var(disp))
D = model[0]/(6*dt) #diffusion coefficient
D_error = model[0]*(1-r2[0])*100/(6*dt)
#print 'diff a', model[0], 'diff b', model[1], 'resid', resid[0], 'r2', r2[0]

time = np.linspace(0,entries,entries)

def plot(x, y):
 plt.figure()
 plt.ylabel(r'y', fontdict=font)
 plt.xlabel(r'x', fontdict=font)
 plt.plot(x, y, 'b', label=r"line")
 plt.legend()
 plt.show()

def plot_energy():
 plt.figure()
 plt.ylabel(r'Energy $(\epsilon)$', fontdict=font)
 plt.xlabel(r'Time ($\Delta t$)', fontdict=font)
 plt.plot(U, 'r', label='Potential energy')
 plt.plot(K, 'g', label='Kinetic energy')
 plt.plot(K+U, 'b', label='Total energy')
 plt.legend(loc=7)
 plt.show()

def plot_pressure():
 plt.figure()
 plt.title(r'Pressure', fontdict=font)
 plt.ylabel(r'$\beta P/ \rho$', fontdict=font)
 plt.xlabel(r'Time ($\Delta t$)', fontdict=font)
 plt.plot(time, P, 'r', label='Pressure', zorder=1)
 plt.xlim(0,len(time))
 plt.hlines(P_avg, 0, len(time), linestyles='dashed', zorder=2)
 plt.legend()
 plt.show()

def plot_pcf():
 histogram = np.loadtxt('histogram.dat')
 bin_size = r_verlet / N_bin
 pcf = np.ones(int(N_bin)-1)
 pcf_x = np.ones(int(N_bin)-1)
 for i in range(1, len(pcf)+1):
  pcf_x[i-1] = bin_size*i
  pcf[i-1] = 2.0*histogram[i-1] / (4*np.pi*bin_size*N*rho*(i*bin_size)**2)
 plt.figure()
 plt.xlim(0,r_verlet)
 plt.ylabel(r'$g(r)$', fontdict=font)
 plt.xlabel(r'$r/\sigma$', fontdict=font)
 plt.hlines(1, 0, r_verlet, linestyles='dashed')
 plt.plot(pcf_x ,pcf, 'r', label=r'Pair correlation function $g(r)$')
 plt.legend()
 plt.show()

def estimated_autocorrelation(x):
 """
 http://stackoverflow.com/q/14297012/190597
 http://en.wikipedia.org/wiki/Autocorrelation#Estimation
 """
 n = len(x)
 variance = x.var()
 x = x-x.mean()
 r = np.correlate(x, x, mode = 'full')[-n:]
 assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
 result = r/(variance*(np.arange(n, 0, -1)))
 return result

def find_nearest(array,value):
 idx = (np.abs(array-value)).argmin()
 return idx, array[idx]


def plot_autocorrelation():
 a = estimated_autocorrelation(U)
 mx = 40
 time = np.linspace(0,mx,mx)
 plt.figure()
 plt.text(mx/2.0-0.05*mx, 0.85, r'$\hat{R}(k)=\frac{1}{(n-k) \sigma^2} \sum_{t=1}^{n-k} (X_t-\mu)(X_{t+k}-\mu) $', fontdict=font)
 plt.ylabel(r'Autocorrelation function', fontdict=font)
 plt.xlabel(r'Time ($\Delta t$)', fontdict=font)
 plt.plot(time, a[0:mx], 'r')
 idx, value = find_nearest(a[0:mx],0.36787944117)
 print idx, value
 plt.annotate(
    '', xy=(0, value), xycoords = 'data',
    xytext = (time[idx], value), textcoords = 'data',
    arrowprops = {'arrowstyle':'<->'})

 plt.annotate(
    r'$\tau$', xy=(time[idx]/2.0, 0.40), xycoords = 'data',
    xytext = (0, 0), textcoords = 'offset points', size=16)
 plt.show()

def plot_temperature():
 plt.xlabel(r'Time ($\Delta t$)', fontdict=font)
 plt.ylabel(r'Temperature ($\epsilon / k$)', fontdict=font)
 plt.plot(T, 'r', label=r'Temperature')
 plt.legend()
 plt.show()

print 'cv', cv_avg, cv_error
print 'mu', mu_avg, mu_error
print 'P ', P_avg, P_error
print 'T ', T_avg, T_error
print 'D ', D, D_error


#Use the following commandos in the terminal to plot
#plot_pressure()
#plot_pcf()
#plot_energy()
#plot_autocorrelation()
#plot(time, K)

