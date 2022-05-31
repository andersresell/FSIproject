import numpy
import numpy as np
import math
import matplotlib.pyplot as plt
from plotting_utilities import *

data = np.loadtxt("../../experiments/Data/Steel_P10_01Load (copy).txt")
t_load = data[:,0]
p_load = data[:,1]
plt.plot(t_load,p_load)
plt.xlabel("t_load")
plt.ylabel("p_load")
#plt.show()
#initial values
#PC
R = 287
p_l = 10.051e5
T_l = 20.581055 + 273.15
p_r = 1005.981628e2
p_amb_pc = p_r
T_r = 20.929594 + 273.15
rho_l = p_l/(R*T_l)
rho_r = p_r/(R*T_r)
print("Pc:")
print("rho_l",rho_l,"p_l",p_l)
print("rho_r",rho_r,"p_r",p_r)
#Al
R = 287
p_l = 9.987e5
T_l = 20.9839 + 273.15
p_r = 100646.68
p_amb_al = p_r
T_r = 21.3694 + 273.15
rho_l = p_l/(R*T_l)
rho_r = p_r/(R*T_r)
print("Al:")
print("rho_l",rho_l,"p_l",p_l)
print("rho_r",rho_r,"p_r",p_r)
#St
R = 287
p_l = 9.928e5
T_l = 21.093750 + 273.15
p_r = 1006.366161e2
p_amb_st = p_r
T_r = 21.449934 + 273.15
rho_l = p_l/(R*T_l)
rho_r = p_r/(R*T_r)
print("St:")
print("rho_l",rho_l,"p_l",p_l)
print("rho_r",rho_r,"p_r",p_r)


output_folder = "output_experiment_free_plate_pc_tube"
p = Plotter(output_folder)
plt.figure()
data = genfromtxt("../../experiments/Data/FreePlate/FP_Al_P10_S01.csv",comments = "#", delimiter=',')
p_S01 = data[:,1]*1e5+4 + p_amb_pc
t_S01 = data[:,0]
ind = np.where(np.logical_and(t_S01>-2,t_S01<3.86))
plt.plot(t_S01,p_S01,'r')
x_S01 = 16.725
time_offset_ms = 29.0831
time,pressure = p.time_history_1D(x_S01,t_a=0.026,t_b=0.034,stride=100)
plt.plot(time*1000-time_offset_ms,pressure,'k')
plt.xlim(-2,3.86)
plt.xlabel(r"""$t\;[\textrm{ms}]$""")
plt.ylabel(r"""$p\;[\textrm{Pa}]$""")
plt.legend(['Experiment','Simulation'])
ax = plt.gca()
ax.locator_params(axis="x", nbins=7)
ax.locator_params(axis="y", nbins=5)
plt.tight_layout()
plt.show()






time_offset_ms = 0.51701
time_offset_ms = 0

p = Plotter("output_experiment_free_plate_pc_box")
t,a,p_wall,u_wall,c_wall = p.piston_fsi_extract_data(t_a=0,t_b=0.00342,stride=10)
ind = np.where(u_wall>-0.5)
plt.figure(2)
plt.plot(t[ind]*1000-time_offset_ms,(a[ind]-a[0])*1000,'k')
plt.figure(3)
plt.plot(t[ind]*1000-time_offset_ms,u_wall[ind],'k')

p = Plotter("output_experiment_free_plate_al_box")
t,a,p_wall,u_wall,c_wall = p.piston_fsi_extract_data(t_a=0,t_b=0.00342,stride=10)
ind = np.where(u_wall>-0.5)
plt.figure(2)
plt.plot(t[ind]*1000-time_offset_ms,(a[ind]-a[0])*1000,'-.k')
plt.figure(3)
plt.plot(t[ind]*1000-time_offset_ms,u_wall[ind],'-.k')

p = Plotter("output_experiment_free_plate_st_box")
t,a,p_wall,u_wall,c_wall = p.piston_fsi_extract_data(t_a=0,stride=10)
ind = np.where(u_wall>-0.5)
plt.figure(2)
plt.plot(t[ind]*1000-time_offset_ms,(a[ind]-a[0])*1000,'--k')
plt.figure(3)
plt.plot(t[ind]*1000-time_offset_ms,u_wall[ind],'--k')

data = genfromtxt("../../experiments/Data/FreePlate/FP_PC_P10_Disp.csv",comments = "#", delimiter=',')
plt.figure(2)
plt.plot(data[:,0],data[:,1],'r')
data = genfromtxt("../../experiments/Data/FreePlate/FP_Al_P10_Disp.csv",comments = "#", delimiter=',')
plt.figure(2)
plt.plot(data[:,0],data[:,1],'r')
data = genfromtxt("../../experiments/Data/FreePlate/FP_St_P10_Disp.csv",comments = "#", delimiter=',')
plt.figure(2)
plt.plot(data[:,0],data[:,1],'r')
plt.xlabel(r"""$t\;[\textrm{ms}]$""")
plt.ylabel(r"""Displacement [mm]""")
#plt.legend(['Sim. PC','Sim. Al','Sim. St','Exp.'])
plt.xlim(0,6)
plt.ylim(0,50)
plt.grid()
plt.locator_params(axis="x", nbins=6)
plt.locator_params(axis="y", nbins=6)
plt.tight_layout()

data = genfromtxt("../../experiments/Data/FreePlate/FP_PC_P10_Vel.csv",comments = "#", delimiter=',')
plt.figure(3)
plt.plot(data[:,0],data[:,1],'r')
data = genfromtxt("../../experiments/Data/FreePlate/FP_Al_P10_Vel.csv",comments = "#", delimiter=',')
plt.figure(3)
plt.plot(data[:,0],data[:,1],'r')
data = genfromtxt("../../experiments/Data/FreePlate/FP_St_P10_Vel.csv",comments = "#", delimiter=',')
plt.figure(3)
plt.plot(data[:,0],data[:,1],'r')
plt.xlabel(r"""$t\;[\textrm{ms}]$""")
plt.ylabel(r"""Velocity [m/s]""")
plt.legend(['Sim. PC','Sim. Al','Sim. St','Exp.'])
#plt.xlim(0,6)
plt.ylim(0,70)
plt.grid()
plt.locator_params(axis="x", nbins=6)
plt.locator_params(axis="y", nbins=6)
plt.tight_layout()


plt.show()