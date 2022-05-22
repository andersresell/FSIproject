import matplotlib.pyplot as plt
import numpy as np

from plotting_utilities import *


#Read header
output_folder = "output_circle_convergence"


p = Plotter(output_folder)
#p.debug_points(232)
#p.debug_animation()
#p.plot_steady_state("p")
#p.contour_animate("p")
p.contour_plot(datatype="p",t0=0.011)
p.plot_convergence()
#plt.show()

error = np.array([1.4485677083333344, 0.8744383169934549, 0.513919375631514, 0.23310947772902227, 0.18219619560158054])
ni = np.array([50, 100, 200, 400, 800])


plt.rcParams.update({
    "font.size": 18,
    #"axes.facecolor": "xkcd:mint green"
    #"axes.facecolor": "lightsteelblue"
    #"lines.linewidth": 2,'
    "axes.labelsize": 22,
    "legend.fontsize": 21,


    #"figure.figsize": (6,6)

})

plt.figure()
plt.plot(ni,error,'*-k')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('ni,nj')
plt.ylabel('Percentage Mass Error')
plt.grid()
plt.tight_layout()
plt.show()