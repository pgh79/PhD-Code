import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
t = np.linspace(0,20,1001)

def x(t,a):
    return np.exp(-a*t)*(np.exp((a-1)*t)-1)/(a-1)

a = np.array([ 0.1,0.2, 0.3, 0.4])

X = np.zeros(shape = (t.size, a.size))

for i in range(a.size):

    X[:,i] = x(t,a[i])


fig, ax = plt.subplots(figsize = (5,3))


ax.plot(t,X)

handles = [r'$k/k_a = {{{A}}}$'.format(A = round(j,2)) for j in a]

legend = ax.legend(handles,
                    framealpha = 0,
                    fontsize = 12
                     )


ax.set_xlabel(r'$k_a t$', fontsize = 12)
ax.set_ylabel(r'$\frac{V}{D} \cdot C(k_a t)$', fontsize = 12)
ax.set_xticks([0,5,10,15,20])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()

plt.savefig('pkcureves.png', dpi = 800)
