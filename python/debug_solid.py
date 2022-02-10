import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
from riemann_exact import riemann_exact
from plotting_utilities import read_header

data = genfromtxt("output_folders/solid_debug/out.csv",comments = "#", delimiter=',')
cell_type = data[:,0]
x = data[:,1]
y = data[:,2]

data = genfromtxt("output_folders/solid_debug_boundary/boundary.csv",comments = "#", delimiter=',')
x_b = data[:,0]
y_b = data[:,1]
plt.plot(x_b,y_b,'.-k')
plt.plot(np.array([x_b[-1],x_b[0]]),np.array([y_b[-1],y_b[0]]),'.-k')

for i in range(0,cell_type.shape[0]):
    if cell_type[i] == 0:
        #fluid
        plt.plot(x[i],y[i],'ob')
    elif cell_type[i] == 1:
        #ghost
        plt.plot(x[i],y[i],'om')
    elif cell_type[i] == 2:
        #solid
        plt.plot(x[i],y[i],'or')

data = genfromtxt("output_folders/solid_debug_intercepts/out.csv",comments = "#", delimiter=',')
x_i = data[:,0]
y_i = data[:,1]
plt.plot(x_i,y_i,'y*')
plt.axis('equal')
plt.show()
