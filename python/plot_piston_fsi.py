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
p.animate_piston_fsi("p",0,10e5)
#p.animate("p")
#p.plot_convergence()
plt.show()