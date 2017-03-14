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


#ax = fig.add_subplot(111, projection='3d')
#ax.plot(data[:,0], data[:,1], data[:,2])

## Calculation of the reward
#reward = np.zeros(data.shape)

dt=0.001;
mass = 1.35
dEv = mass/dt*data[:-1,3]*(data[1:,3]-data[:-1,3])
dEp = mass/dt*(data[1:,2] - data[:-1,2])

dEtot = dEv + dEp

plt.figure()

lin_energytot, = plt.plot(dEtot, label='Energie totale')
lin_energypot, = plt.plot(dEp, label='Energie potentielle')
lin_energycin, = plt.plot(dEv, label='Energie Cine')

#plt.legend()handles=[lin_energytot, lin_energypot,lin_energycin])

plt.show()

