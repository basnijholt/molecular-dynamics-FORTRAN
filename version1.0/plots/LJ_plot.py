import numpy as np
import matplotlib.pyplot as plt


font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

x = np.linspace(0.0, 4.0, 200)
y = 4*(1/x**12-1/x**6)
plt.ylim(-1.5,1.5)
plt.xlim(0,4)
plt.plot(x, y, 'k')
plt.title('Lennard-Jones potential', fontdict=font)
plt.text(2, 0.65, r'$\phi_{LJ}(r) = 4 \left[  \left(  \frac{1}{r}\right)^{12} - \left(  \frac{1}{r}\right)^{6} \right]$', fontdict=font)
plt.xlabel(r'r ($\sigma$)', fontdict=font)
plt.ylabel(r'Potential ($\epsilon$)', fontdict=font)


plt.hlines(0, 1, 4, linestyles='dashed')
plt.vlines(2**(1.0/6.0), 0, -1, linestyles='dashed')
plt.hlines(-1, 0, 4, linestyles='dashed')

text_grote=16
# de pijl voor sigma
plt.annotate(
    '', xy=(1, 0), xycoords = 'data',
    xytext = (0, 0), textcoords = 'data',
    arrowprops = {'arrowstyle':'<->'})

plt.annotate(
    r'$\sigma$', xy=(0.5, 0.05), xycoords = 'data',
    xytext = (0, 0), textcoords = 'offset points', size=text_grote)


# de pijl voor epsilon
plt.annotate(
    '', xy=(0.5, 0), xycoords = 'data',
    xytext = (0.5, -1), textcoords = 'data',
    arrowprops = {'arrowstyle':'<->'})

plt.annotate(
    r'$\epsilon$', xy=(0.55, -0.5), xycoords = 'data',
    xytext = (0, 0), textcoords = 'offset points', size=text_grote)

# r_cut off
plt.annotate(
    '', xy=(3, -0.5), xycoords = 'data',
    xytext = (2.5, -0.01), textcoords = 'data',
    arrowprops = {'arrowstyle':'<-'})

plt.annotate(
    r'$r_{cut}$', xy=(3.05, -0.55), xycoords = 'data',
    xytext = (0, 0), textcoords = 'offset points', size=text_grote)

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()
