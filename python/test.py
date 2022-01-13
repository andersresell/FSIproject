
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt


n_timesteps = 1500

timestride = 10

for n in range(0,n_timesteps+1):
    if n % timestride == 0:
        my_data = genfromtxt("fvm_output/fvm_out_t"+str(n)+".csv",comments = "#", delimiter=',')
        if n == 0:
            ni = int(my_data[0,0])
            nj = int(my_data[0,1])
            L_x = my_data[0,2]
            L_y = my_data[0,3]
            dx = L_x/ni
            dy = L_y/nj
            x = np.linspace(dx/2,L_x-dx/2,nj)
            y = np.linspace(dy/2,L_y-dy/2,nj)
        my_data = my_data[1::,:]
        E = my_data[:,3]
        E = np.transpose(E.reshape((ni,nj)))

        plt.clf()
        cs = plt.contourf(x,y,E)
        plt.pause(0.000001)
    print(n)
plt.show()
