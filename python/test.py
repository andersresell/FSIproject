
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
            #x_stride = np.arange(0,ni,int(ni/100))
            #y_stride = np.arange(0,nj,int(nj/100))
            x = np.linspace(dx/2,L_x-dx/2,ni)
            y = np.linspace(dy/2,L_y-dy/2,nj)
            levels = np.linspace(0,0.5e6,100)
        my_data = my_data[1::,:]
        p = my_data[:,3]
        p = np.transpose(p.reshape((ni,nj)))

        plt.clf()
        cs = plt.contourf(x,y,p, levels=levels, cmap=plt.get_cmap('hot'))
        plt.axis('equal')
        plt.pause(0.00001)
    print(n)
plt.show()
