import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D

# Lecture du fichier :
fname = "DATA/energy.txt"
data = pd.read_csv(fname,sep = ' ')

plt.figure()

lin_dEc = plt.plot(data["dEc"], label='dEc')
lin_dEp = plt.plot(data["dEp"], label='dEp')
lin_dEtot = plt.plot(data["dEtot"], label='dEtot')

plt.legend()

plt.show()

