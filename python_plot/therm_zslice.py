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
ax = fig.gca(projection='3d')
#fig.suptitle('Allen Model ', fontsize=20)
plt.xlabel('Distance along X direction[m]', fontsize=10)
plt.ylabel('Distance along Y direction[m]', fontsize=10)
plt.title('Updraft velocity[m/s]', fontsize=15)
cs=ax.plot_trisurf(data_t["x"], data_t["y"], data_t["updraft"], cmap=cm.jet, linewidth=0.2)
#proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) 
#    for pc in cs.collections]

#plt.legend(["range(2-3)", "range(3-4)", "range(4-6)"])
#plt.show()
#fig.savefig('test.jpg')
plt.show()	
