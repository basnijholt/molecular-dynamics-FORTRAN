import numpy as np
import matplotlib.pyplot as plt


font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }


N = 4
L = 5

x = np.zeros(N)
y = np.zeros(N)

L_min = L/10.0
L_max = L*9.0/10.0
L_box = L_max - L_min


for i in range(0,N):
 x[i] = np.random.uniform(L_min,L_max)
 y[i] = np.random.uniform(L_min,L_max)  


plt.subplot(2, 1, 1)
plt.ylim(-1,6)
plt.xlim(-1,6)
plt.plot(x, y, 'ro', label='Particles')
plt.plot(-5, -5, 'bo', label='Mirrored Particles')
plt.title('2D box with periodic boundary conditions', fontdict=font)
van = -1
tot = 6
plt.hlines(L_min, van, tot, linestyles='dashed')
plt.hlines(L_max, van, tot, linestyles='dashed')
plt.vlines(L_min, van, tot, linestyles='dashed')
plt.vlines(L_max, van, tot, linestyles='dashed')
for i in range(0,N-1):
 for j in range(i+1,N):
  r = np.sqrt( (x[i]-x[j])**2 + (y[i]-y[j])**2 )
  bc_xj = 0
  bc_yj = 0
  bc_xi = 0
  bc_yi = 0
  if r < 0.5*L_box:
   plt.annotate('', xy=(x[i],y[i]), xycoords='data', xytext=(x[j],y[j]), textcoords='data', arrowprops = {'arrowstyle':'<->'})
  elif np.abs(x[i]-x[j]) > 0.5*L_box: #afstand groter dan half
   if x[i] < 0.5*L:
    bc_xj = -1
   else:
    bc_xj = 1
   if x[j] < 0.5*L:
    bc_xi = -1
   else:
    bc_xi = 1
  elif np.abs(y[i]-y[j]) > 0.5*L_box: #afstand groter dan half
   if y[i] < 0.5*L:
    bc_yj = -1
   else:
    bc_yj = 1
   if y[j] < 0.5*L:
    bc_yi = -1
   else:
    bc_yi = 1
  xj = x[j] + bc_xj * L_box
  yj = y[j] + bc_yj * L_box
  xi = x[i] + bc_xi * L_box
  yi = y[i] + bc_yi * L_box
  plt.annotate('', xy=(xj,yj), xycoords='data', xytext=(x[i],y[i]), textcoords='data', arrowprops = {'arrowstyle':'<->'})
  plt.annotate('', xy=(xi,yi), xycoords='data', xytext=(x[j],y[j]), textcoords='data', arrowprops = {'arrowstyle':'<->'})
  if ((xi > L_max) | (xi < L_min)) | ((yi > L_max) | (yi < L_min)):
   plt.plot(xi, yi, 'bo')
  if ((xj > L_max) | (xj < L_min)) | ((yj > L_max) | (yj < L_min)):
   plt.plot(xj, yj, 'bo')

plt.legend()

plt.subplot(2, 1, 2)
plt.ylim(-1,6)
plt.xlim(-1,6)
plt.plot(x, y, 'ro', label='Particles')
plt.title('2D box with hard boundary conditions', fontdict=font)
plt.xlabel(r'x-axis ($\sigma$)', fontdict=font)
plt.ylabel(r'y-axis ($\sigma$)', fontdict=font)
plt.hlines(L_min, L_min, L_max, linestyles='-')
plt.hlines(L_max, L_min, L_max, linestyles='-')
plt.vlines(L_min, L_min, L_max, linestyles='-')
plt.vlines(L_max, L_min, L_max, linestyles='-')
for i in range(0,N-1):
 for j in range(i+1,N):
  plt.annotate('', xy=(x[i],y[i]), xycoords='data', xytext=(x[j],y[j]), textcoords='data', arrowprops = {'arrowstyle':'<->'})

plt.show()
