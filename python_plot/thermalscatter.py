import sys, os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

filename="DATA/wind.txt"
data = pd.read_csv(filename,sep = ' ')
data_t = data[data["t"]==400]
#data_t_z = data_t[data_t["z"]==900]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

sp = ax.scatter(data_t["x"], data_t["y"], data_t["z"], s=20, c=data_t["updraft"])
plt.colorbar(sp)

plt.show()
