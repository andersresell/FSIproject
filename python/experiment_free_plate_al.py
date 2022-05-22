import numpy
import numpy as np
import math
import matplotlib.pyplot as plt
from plotting_utilities import *

output_folder = "output_experiment_free_plate_al_box"
p = Plotter(output_folder)
#p.debug_points(0)
#plt.show()
p.animate_piston_fsi(datatype="p",bottom_level=0,top_level=10e5)
plt.show()

#data = genfromtxt("../../experiments/Data/FreePlate/FP_Al_P10_Force.csv",comments = "#", delimiter=',')
#plt.figure()
#plt.plot(data[:,0],data[:,1])

p_amb = 100646.68

output_folder = "output_experiment_free_plate_al_tube"
p = Plotter(output_folder)
#p.animate_1D(datatype="p",bottom_level=0,top_level=1e6)

plt.figure()
data = genfromtxt("../../experiments/Data/FreePlate/FP_Al_P10_S01.csv",comments = "#", delimiter=',')
p_S01 = data[:,1]*1e5+4 + p_amb
t_S01 = data[:,0]
ind = np.where(np.logical_and(t_S01>-2,t_S01<3.86))
plt.plot(t_S01,p_S01,'r')
x_S01 = 16.725
time_offset_ms = 28.9362
time,pressure = p.time_history_1D(x_S01,t_a=0.026,t_b=0.034)
plt.plot(time*1000-time_offset_ms,pressure,'k')
plt.xlim(-2,3.86)
plt.xlabel(r"""$t\;[\textrm{ms}]$""")
plt.ylabel(r"""$p\;[\textrm{Pa}]$""")
plt.legend(['Experiment','Simulation'])
plt.locator_params(axis="y", nbins=5)
plt.tight_layout()
#plt.show()
#plt.figure()
#p.plot_1D(datatype="p",bottom_level=0,top_level=1e6,t_goal=0.0275)
#plt.show()
output_folder = "output_experiment_free_plate_al_box"
p = Plotter(output_folder)
p.animate_piston_fsi(datatype="p",bottom_level=0,top_level=20e5)
plt.figure()
p.debug_points(0)
plt.show()
#x_S01 = 0.1835
#time,pressure = p.time_history_1D(x_S01)#,t_b=0.001)
#plt.figure()
#plt.plot(time*1000,pressure,'k')


t,a,p_wall,u_wall = p.piston_fsi_extract_data(t_a=0,t_b=0.00338)
data = genfromtxt("../../experiments/Data/FreePlate/FP_Al_P10_Disp.csv",comments = "#", delimiter=',')
plt.figure()
plt.plot(data[:,0],data[:,1],'r')
time_offset_ms = 0.51701
plt.plot(t*1000,(a-a[0])*1000,'k')
plt.xlabel(r"""$t\;[\textrm{ms}]$""")
plt.ylabel(r"""Displacement [mm]""")
plt.legend(['Experiment','Simulation'])
plt.tight_layout()

data = genfromtxt("../../experiments/Data/FreePlate/FP_Al_P10_Vel.csv",comments = "#", delimiter=',')
plt.figure()
plt.plot(data[:,0],data[:,1])
plt.figure()
plt.plot(t*1000,u_wall)



plt.show()