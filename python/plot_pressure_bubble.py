import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import genfromtxt
from analytical_solutions import riemann_exact
from plotting_utilities import *


#Read header

output_folder = "output_pressure_bubble"

p = Plotter(output_folder)
#p.debug_points()
#p.plot_steady_state("p")
p.contour_animate("p",1e4,5e5)
plt.show()