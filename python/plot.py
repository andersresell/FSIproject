
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt


fvm_output_folder = "fvm_output"

#Read header
my_data = genfromtxt(fvm_output_folder+"/header.csv",comments = "#", delimiter=',')
ni = int(my_data[0])
nj = int(my_data[1])
L_x = my_data[2]
L_y = my_data[3]
n_timesteps = int(my_data[4])
write_stride = int(my_data[5])
dx = L_x/ni
dy = L_y/nj
x = np.linspace(dx/2,L_x-dx/2,ni)
y = np.linspace(dy/2,L_y-dy/2,nj)

levels = np.linspace(0,0.5e6,100)

#animate
for n in range(0,n_timesteps+1):
    if n % write_stride == 0:
        my_data = genfromtxt(fvm_output_folder+"/fvm_out_t"+str(n)+".csv",comments = "#", delimiter=',')
        p = my_data[:,3]
        p = np.transpose(p.reshape((ni,nj)))

        plt.clf()
        cs = plt.contourf(x,y,p, levels=levels, cmap=plt.get_cmap('hot'))
        plt.axis('equal')
        plt.pause(0.00001)


    print(n)

plt.show()
