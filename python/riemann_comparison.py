

import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
from riemann_exact import riemann_exact
from plotting_utilities import read_header

#Read header
fvm_output_folder = "output_riemann_sod"
ni,nj,L_x,L_y,write_stride,n_timesteps, t_end, dx,dy, x,y = read_header(fvm_output_folder)


#Plotting the last step
levels = np.linspace(0,3,100)
my_data = genfromtxt("output_folders/"+fvm_output_folder+"/fvm_out_t"+str(n_timesteps)+".csv",comments = "#", delimiter=',')
rho = my_data[:,0]
rho = np.transpose(rho.reshape((ni,nj)))
u = my_data[:,1]
u = np.transpose(u.reshape((ni,nj)))
p = my_data[:,3]
p = np.transpose(p.reshape((ni,nj)))


figure, axis = plt.subplots(2,2)
axis[0,0].plot(x,rho[int(nj/2),:],'.')
axis[1,0].plot(x,u[int(nj/2),:],'.')
axis[0,1].plot(x,p[int(nj/2),:],'.')
#plt.plot(x,p[int(nj/2),:],'.') #plotting the pressure along y = L_y/2

#Sod's problem
rho_l = 3
u_l = 0
p_l = 3
rho_r = 1
u_r = 0
p_r = 1

#Stationary contact discontinuity problem
#rho_l = 3
#u_l = 1
#p_l = 1
#rho_r = 1
#u_r = 1
#p_r = 1

rho,u,p = riemann_exact(rho_l,u_l,p_l,rho_r,u_r,p_r,ni,L_x,t_end)
#plt.plot(x,p)

axis[0,0].plot(x,rho)
axis[0,0].set_xlabel("x")
axis[0,0].set_ylabel("Density")

axis[1,0].plot(x,u)
axis[1,0].set_xlabel("x")
axis[1,0].set_ylabel("Velocity")

axis[0,1].plot(x,p)
axis[0,1].set_xlabel("x")
axis[0,1].set_ylabel("Pressure")

fig = plt.gcf()
fig.set_size_inches(8, 6)

plt.show()
