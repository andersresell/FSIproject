import numpy
import numpy as np
import math
import matplotlib.pyplot as plt
from plotting_utilities import *

R = 287
T_a = 22.8+273.15
p_a = 101858
rho_a = p_a/(R*T_a)
T_d = 22.3+273.15
p_d = 10e5 + p_a
#p_d = 2.3e5+p_a
rho_d = p_d/(R*T_d)
print("T_a",T_a)
print("p_a",p_a)
print("rho_a",rho_a)
print("T_d",T_d)
print("p_d",p_d)
print("rho_d",rho_d)

file = "../../experiments/Square/S_03/S_03.txt"
data = numpy.loadtxt(file,skiprows=19)
t = data[:,0]*1e-3
P08_2= data[:,8]
P08_1 = data[:,7]
P08_2 = P08_2 - np.mean(P08_2[0:20])
P08_1 = P08_1 - np.mean(P08_1[0:20])
P08_2 = P08_2*1e5 + p_a
P08_1 = P08_1*1e5 + p_a
x_P08_2 = 4.1
x_P08_1 = 4.2
y_center = 0.29/2

plt.plot(t,P08_2)
plt.plot(t,P08_1)


#file = "../../experiments/Square/S_03/S_03_synced.csv"
#data = genfromtxt(file, delimiter=',',skip_header=2)
#P08_2 = data[:,5]*1e5 + p_a
#t = np.zeros(P08_2.shape)
#for i in range(0,t.size):
#    t[i] = i/37000
#plt.plot(t,P08_2)


output_folder = "output_experiment_square_3_tube"
p = Plotter(output_folder)
time,pressure = p.time_history(x_P08_2,y_center)
plt.plot(time,pressure,'.')
#plt.show()
time,pressure = p.time_history(x_P08_1,y_center)
plt.plot(time,pressure)
plt.xlabel("t[s]")
plt.ylabel("p[Pa]")
plt.tight_layout()
#p.animate_1D("p",0,12e5)
#p.animate_1D("u",0,300)
#p.animate_1D("rho",0,10)
#p.contour_animate("p")#,100000,150000)
plt.show()

#p.contour_animate("p")

output_folder = "output_experiment_square_3_box"
p = Plotter(output_folder)

#p.debug_points_last()


#p.animate_1D("p")#,0,400000)
#p.contour_animate("p")#,100000,150000)
#p.animate_1D("p",0,350000)
#p.schlieren_last()
#p.schlieren_animate()
#p.schlieren_plot_all()
#p.schlieren_plot_time_interval(t0=0.000235,show_timestep=False)



plt.show()