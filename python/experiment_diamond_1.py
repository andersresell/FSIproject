import numpy
import numpy as np
import math
import matplotlib.pyplot as plt
from plotting_utilities import *

R = 287
T_a = 23.03+273.15
p_a = 101890
rho_a = p_a/(R*T_a)
T_d = 21.35+273.15
p_d = 2.16e5 + p_a
#p_d = 2.3e5+p_a
rho_d = p_d/(R*T_d)
dt = 0.002e-3
print("T_a",T_a)
print("p_a",p_a)
print("rho_a",rho_a)
print("T_d",T_d)
print("p_d",p_d)
print("rho_d",rho_d)

t_P08_2_sim = 0.008724*1e3
t_P08_2_exp = 0.049999*1e3
delta_t_P08_2 = t_P08_2_exp-t_P08_2_sim

t_P08_1_sim = 0.008950294
t_P08_1_exp = 0.05023584
delta_t_P08_1 = t_P08_1_exp-t_P08_1_sim




file = "../../experiments/Diamond/D_01/D_01.txt"
data = numpy.loadtxt(file,skiprows=19)
t = data[:,0]#*1e-3
P08_2= data[:,8]
P08_1 = data[:,7]
P08_2 = P08_2 - np.mean(P08_2[0:20])
P08_1 = P08_1 - np.mean(P08_1[0:20])
P08_2 = P08_2*1e5 + p_a
P08_1 = P08_1*1e5 + p_a
x_P08_2 = 4.1
x_P08_1 = 4.2
y_center = 0.29/2

plt.plot(t-delta_t_P08_2,P08_2,'r')
#plt.plot(t-delta_t_P08_1,P08_1)



output_folder = "output_experiment_diamond_1_tube"
p = Plotter(output_folder)
time,pressure = p.time_history_1D(x_P08_2)
plt.plot(time*1e3,pressure,'k')
#time,pressure = p.time_history_1D(x_P08_1)
#plt.plot(time,pressure)
#x_test = 5.075#center of test section
#time,pressure = p.time_history_1D(x_test)
#plt.plot(time,pressure)

file = "../../experiments/Diamond/D_01/D_01_synced.csv"
data = genfromtxt(file, delimiter=',',skip_header=2)
P08_2 = data[:,5]*1e5 + p_a
t = np.zeros(P08_2.shape)
for i in range(0,t.size):
    t[i] = i/37000
#t = t*1e3
#plt.plot(t,P08_2)
plt.xlim(8,12)
plt.xlabel("t [ms]")
plt.ylabel("p [Pa]")
plt.legend(["Experiment","Simulation"])
ax = plt.gca()
ax.locator_params(axis="x", nbins=5)
ax.locator_params(axis="y", nbins=5)
#plt.ylim(90000,190000)
plt.tight_layout()
#plt.show()
plt.figure()
p.plot_1D("p",t_goal=t_P08_2_sim*1e-3,show_timestep=False)
#p.animate_1D("p",0,320000)
#p.contour_animate("p")#,100000,150000)
#plt.show()


#p.contour_animate("p")

output_folder = "output_experiment_diamond_1_box"
p = Plotter(output_folder)




#p.animate_1D("p")#,0,400000)
#p.contour_animate("p")#,100000,150000)
#p.animate_1D("p")
#p.schlieren_last()
#p.schlieren_animate()
#p.schlieren_plot_all()
#p.schlieren_plot_time_interval()



plt.show()