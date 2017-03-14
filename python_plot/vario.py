#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:49:55 2017

@author: adrien
"""


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

data_coordonnees = data[:,[0,1,2]]

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot(data[:,0], data[:,1], data[:,2])

## Calculation of the reward
#reward = np.zeros(data.shape)

mass = 1.35
energy_cinetique = 1/2 * mass*(data[1:,3]*data[1:,3] - data[:-1,3]*data[:-1,3])
energy_potentiel = 9.81 *mass* (data[1:,2] - data[:-1,2])

print(energy_cinetique.shape)
print(energy_potentiel.shape)

energy_tot = energy_cinetique + energy_potentiel

ax2 = fig.add_subplot(111)
ax2.plot(energy_tot)
ax2.plot(energy_potentiel)
ax2.plot(energy_cinetique)
#ax2.plot(1/2 * mass*(data[1:200,3]*data[1:200,3] - data[2:201,3]*data[2:201,3]))

print(energy_cinetique)

plt.show()

