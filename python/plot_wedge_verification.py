import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import genfromtxt
from analytical_solutions import riemann_exact
from plotting_utilities import *


#Read header
output_folder = "wedge_verification_attached"
#output_folder = "wedge_verification_detached"

p = Plotter(output_folder)
#p.contour_animate("p")
p.contour_plot("M",-1,True)
p.plot_convergence()
print(p.probe("M",9,14))
print(p.probe("M",4,7))
print(p.probe("M",4.67,6.23))
plt.show()


