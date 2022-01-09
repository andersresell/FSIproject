
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt

Lx = 1
Ly = 1
ni = 100
nj = 100
dx = Lx/ni
dy = Ly/nj
n_timesteps = 100
x = np.linspace(dx/2,Lx-dx/2,ni)
y = np.linspace(dy/2,Ly-dy/2,nj)

t = 0
for n in range(0,n_timesteps+1):
    my_data = genfromtxt("fvm_output/fvm_out_t"+str(n)+".csv",delimiter=',',skip_header=1)
    E = my_data[:,3]
    E = E.reshape((ni,nj))
    plt.clf()
    cs = plt.contourf(x,y,E)
    plt.pause(0.000001)
    print(n)
plt.show()
