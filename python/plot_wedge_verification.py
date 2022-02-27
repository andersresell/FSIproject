import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import genfromtxt
from riemann_exact import riemann_exact
from plotting_utilities import *


#Read header
output_folder = "wedge_verification_attached"
#output_folder = "wedge_verification_detached"

p = Plotter(output_folder)
#p.debug_points()
p.plot_steady_state("p")
p.animate("p")

plt.show()


