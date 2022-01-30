
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
from riemann_exact import riemann_exact
from plotting_utilities import read_header

#Read header
fvm_output_folder = "output_riemann_sod"
ni,nj,L_x,L_y,write_stride,n_timesteps, t_end, dx,dy, x,y = read_header(fvm_output_folder)



levels = np.linspace(0,3,100)

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
