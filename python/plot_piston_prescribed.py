import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import genfromtxt
from analytical_solutions import riemann_exact
from plotting_utilities import *


#Read header
output_folder = "output_piston_prescribed"


p = Plotter(output_folder)

p.contour_plot("p")
p.plot_piston_prescribed(500)
p.plot_convergence()
plt.show()

