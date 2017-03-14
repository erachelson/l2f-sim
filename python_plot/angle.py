
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D

"""
Dans cette partie on plot les résultats des fichiers de sorties
"""

"""
on utilise le script comme ca
./data_plot.py nom_du_fichier
avec le nom du fichier un fichier .txt avec lignes et colonnes
"""

## Lecture du fichier :
fname = "DATA/data_plane.txt"
data = np.loadtxt(fname,dtype=float)



plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot(data[:,0], data[:,1], data[:,2])

## Calculation of the reward
#reward = np.zeros(data.shape)


alpha = data[:,6]
beta = data[:,7]
gamma = data[:,4]
khi = data[:,5]
sigma = data[:,8]

linalp, = plt.plot(alpha, label='alpha')
linbet, = plt.plot(beta, label='beta')
lingam, = plt.plot(gamma, label='gamma')
linsig, = plt.plot(sigma, label='sigma')
linkhi, = plt.plot(khi, label='khi')

plt.legend(handles=[linalp, linbet,lingam,linsig,linkhi])

plt.show()