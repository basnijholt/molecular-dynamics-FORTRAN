import numpy as np
import matplotlib.pyplot as plt
# files I want to read 
file = ['t.dat', 'u.dat', 'p.dat', 'kin_energy.dat']
data_block_size = 50
counter = 0

print 'time_steps'
print np.loadtxt('time_end.dat') 

print 'rho'
print np.loadtxt('rho.dat')

print 't0'
print np.loadtxt('t0.dat')

for i in range(np.size(file)):
 the_array = np.loadtxt(file[i])
 cor_len = np.size(the_array)-np.size(the_array)%data_block_size
 right_size_array = the_array[:cor_len].reshape(data_block_size,-1)
 averages = np.mean(right_size_array, axis = 0)
 average = np.mean(averages, axis = 0)
 var = np.std(averages)/np.sqrt(np.size(averages))
 print file[i]
 print average
 print var

mu = np.loadtxt('mu.dat')
print 'mu'
print np.mean(mu, axis = 0)
print np.std(mu)

# de least-square methode voor de fit van de diffusie coefficient
y = np.loadtxt('disp_sq.dat')
n = np.size(y)
x = np.arange(n)
A = np.vstack([x, np.ones(n)]).T
model, resid = np.linalg.lstsq(A, y)[:2]
r2 = 1 - resid / (np.size(y) * np.var(y))
print 'diff a'
print model[0]
print 'diff b'
print model[1]
print 'resid'
print resid[0]
print 'r2'
print r2[0]

# laad de Cv
cv = np.loadtxt('cv.dat')
print 'cv'
print cv

# laad de eq_time
eq_time = np.loadtxt('eq_time.dat')
print 'eq_time'
print eq_time

print 'delta T in %'


