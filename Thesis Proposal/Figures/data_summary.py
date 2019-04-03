#!/Users/demetri/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import cycler #For color cycle

#For tex fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')



df = pd.read_csv('~/Documents/PhD COde/2019-02-07 LHRD Results/Data/apixiban_regression_data.csv')

grid = sns.FacetGrid(data = df, col='Sex', height = 3, aspect = 1)

color = plt.cm.RdBu(np.linspace(0.1, 0.9,2))
mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)


grid.map_dataframe(sns.lineplot,
                   hue= 'Group',
                   x = 'Time',
                   y = 'Concentration_scaled',
                   ci = 'sd',
                   linewidth = 2
                  )

grid.axes[0,0].set_ylabel('Concentration (mg/L)', fontsize = 12)
grid.axes[0,0].set_xlabel('Time (Hours)',fontsize = 12)
grid.axes[0,1].set_xlabel('Time (Hours)',fontsize = 12)

plt.legend(frameon = False)
plt.tight_layout()
grid.fig.savefig('data_summary.png', dpi = 800, transparent = True)



# plt.tight_layout()
# plt.savefig('data_summary.png', dpi = 800, transparent = True)
