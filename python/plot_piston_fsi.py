import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import genfromtxt
from analytical_solutions import riemann_exact
from plotting_utilities import *


#Read header
output_folder = "output_piston_fsi"


p = Plotter(output_folder)
#p.debug_points(10)

times = np.array([0,20,25,29,33])*1e-3
for i in range(0,times.size):
    p.plot_piston_fsi(datatype="p",t0=times[i],bottom_level=0,top_level=7e5)
#p.contour_plot("p")
#p.contour_animate("p")
#for n in range(1000,1900):
 #   if n%(4*60) == 0:
  #      plt.figure()
   #     p.plot_piston_fsi("p",n,0,10e5)

#p.animate_piston_fsi("p")
#plt.figure()
#p.plot_piston_fsi("p",1903,0,10e5)
#p.contour_animate("p",0,10e5)
#p.contour_plot("p",-1,False)
#p.plot_convergence()
plt.show()