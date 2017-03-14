# -*- coding: utf-8 -*-
"""
    Created on Wed Jan 18 15:38:06 2017
    
    @author: adrien Bufort
    """

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
import sys, os
from matplotlib import cm

"""
    Dans cette partie on plot les résultats des fichiers de sorties
    """

"""
    on utilise le script comme ca
    ./data_plot.py nom_du_fichier
    avec le nom du fichier un fichier .txt avec lignes et colonnes
    """

## Lecture du fichier :
fname = "../data_plane.txt"
data = np.loadtxt(fname,dtype=float)

data_coordonnees = data[:,[0,1,2]]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(data[:,0], data[:,1], data[:,2])

filename="../wind.txt"
data = pd.read_csv(filename,sep = ' ')
data_t = data[data["t"]==400]
data_t_z = data_t[data_t["z"]==410]

#fig.suptitle('Allen Model ', fontsize=20)
plt.xlabel('Distance along X direction[m]', fontsize=10)
plt.ylabel('Distance along Y direction[m]', fontsize=10)
plt.title('Updraft velocity[m/s]', fontsize=15)
cs=ax.plot_trisurf(data_t_z["x"], data_t_z["y"], data_t_z["updraft"], cmap=cm.jet, linewidth=0.2)
#proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0])
#    for pc in cs.collections]

#plt.legend(["range(2-3)", "range(3-4)", "range(4-6)"])
#plt.show()
#fig.savefig('test.jpg')

plt.show()
