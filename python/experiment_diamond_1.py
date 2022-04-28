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
rho_d = p_d/(R*T_d)
dt = 0.002e-3
print("T_a",T_a)
print("p_a",p_a)
print("rho_a",rho_a)
print("T_d",T_d)
print("p_d",p_d)
print("rho_d",rho_d)

file = "../../experiments/Diamond/D_01/D_01.txt"
data = numpy.loadtxt(file,skiprows=19)
t = data[:,0]*1e-3
P10_2= data[:,10]
P10_1 = data[:,9]
P10_2 = P10_2 - np.mean(P10_2[0:20])
P10_1 = P10_1 - np.mean(P10_1[0:20])
P10_2 = P10_2*1e5 + p_a
P10_1 = P10_1*1e5 + p_a
x_P10_2 = 0.26+3.84
x_P10_1 = 0.26+3.94
y_center = 0.29/2



output_folder = "output_experiment_diamond_1"
p = Plotter(output_folder)
p.schlieren_experiment()
#time,pressure = p.time_history(x_P10_2,y_center)
plt.figure()
#plt.plot(time,pressure)
plt.plot(t,P10_2)
plt.xlabel("t[s]")
plt.ylabel("p[Pa]")
plt.tight_layout()

#p.animate_1D("p",0,400000)


#p.contour_animate("p")#,100000,300000)



plt.show()