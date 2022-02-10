import numpy as np
from numpy import genfromtxt

def read_header(fvm_output_folder):
    my_data = genfromtxt("output_folders/"+fvm_output_folder+"/header.csv",comments = "#", delimiter=',')
    #my_data = genfromtxt("fvm_output_folder"/header.csv",comments = "#", delimiter=',')
    ni = int(my_data[0])
    nj = int(my_data[1])
    L_x = my_data[2]
    L_y = my_data[3]
    write_stride = int(my_data[4])
    n_timesteps = int(my_data[5])
    t_end = my_data[6]
    dx = L_x/ni
    dy = L_y/nj
    x = np.linspace(dx/2,L_x-dx/2,ni)
    y = np.linspace(dy/2,L_y-dy/2,nj)
    return ni,nj,L_x,L_y,write_stride,n_timesteps, t_end, dx,dy, x,y