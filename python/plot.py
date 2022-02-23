
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
from riemann_exact import riemann_exact
from plotting_utilities import *

#Read header
#fvm_output_folder = "wedge_verification_attached"
output_folder = "wedge_verification_detached"
ni,nj,L_x,L_y,write_stride,n_timesteps, t_end, dx,dy, x,y, n_static_solids, n_movable_solids = read_headers(output_folder)



#levels = np.linspace(0.25e5,3e5,100)

#animate
for n in range(0,n_timesteps+1):
    if n % write_stride == 0:
        my_data = genfromtxt("output_folders/"+output_folder+"/fvm_output_t"+str(n)+".csv",comments = "#", delimiter=',')
        p = my_data[:,3]
        p = np.transpose(p.reshape((ni,nj)))

        plt.clf()
        cs = plt.contourf(x,y,p, 100, cmap=plt.get_cmap('hot'))
        plot_solids(output_folder,n,n_static_solids,n_movable_solids)
        plt.axis('equal')

        plt.pause(0.00001)


    print(n)
plt.show()
