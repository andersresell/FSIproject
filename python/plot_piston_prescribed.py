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
#p.debug_points(2)
#p.debug_points_last()
#p.debug_animation()
#p.plot_steady_state("p")
#p.animate_1D("p")
#p.animate("p")
p.plot_piston_prescribed(1,100)
p.plot_convergence()
plt.show()