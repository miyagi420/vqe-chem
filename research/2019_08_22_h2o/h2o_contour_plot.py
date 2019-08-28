import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

fig_name = 'contour_plot'

x = np.linspace(0.6,1.5,12)
y = np.linspace(30,180,12)
X, Y = np.meshgrid(x, y)
Z = np.loadtxt('resources/contour_results.txt').reshape(12,12)

cmap= plt.cm.get_cmap("CMRmap", 55)
levels = np.linspace(-75.04,-73,55)

contours = plt.contour(X, Y, Z, 
            levels=levels, 
            colors='white',
            alpha=0.7,
            linewidths=0.5
            )

plt.contourf(X, Y, Z, 
            levels=levels, 
            cmap=cmap, 
            alpha=1
            )

plt.colorbar();
plt.scatter(X,Y,marker='.',s=1,color='orange')

plt.gca().clabel(contours, contours.levels[0:8], inline=True, fontsize=8,
                 fmt="%2.3f", use_clabeltext=True)
plt.title('Geometry-Energy Graph for the H$_{2}$O Molecule')
plt.xlabel('Length of each OH bond (Angstroms)')
plt.ylabel('Angle between OH bonds (degrees)')
plt.savefig('figures/{}.png'.format(fig_name))
