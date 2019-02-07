
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cycler #For color cycle

#For tex fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

t = np.linspace(0,12,1001)

def PKCURVE(t,a):

    #Non dimensionalized solution to 1 compartment model
    return np.exp(-a*t)*(np.exp((a-1)*t)-1)/(a-1)

#Values of the non-dimensional parameter
a = np.array([0.15, 0.2, 0.3, 0.4])

#Initialize a place to store the curves
X = np.zeros(shape = (t.size, a.size))

#Set up the color cycle
n = a.size
color = plt.cm.RdBu(np.linspace(0.1, 0.9,n))
mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

#Populate X with the curves
for i in range(a.size):

    X[:,i] = PKCURVE(t,a[i])


fig, ax = plt.subplots(figsize = (5,3))
ax.plot(t,X)

handles = [r'$k/k_a = {{{A}}}$'.format(A = round(j,2)) for j in a]
legend = ax.legend(handles,
                    framealpha = 0,
                    fontsize = 12,
                    loc = 'upper right'
                     )


ax.set_xlabel(r'$k_a t$', fontsize = 12)
ax.set_ylabel(r'$\frac{V}{D} \cdot C(k_a t)$', fontsize = 12)
ax.set_xticks(np.arange(t.min(),t.max()+1,4))
ax.set_yticks(np.arange(0,1,0.2))

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig('pkcureves.png', dpi = 800, transparent = True)
