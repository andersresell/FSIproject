import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
from analytical_solutions import riemann_exact
from plotting_utilities import *

p = Plotter("output_static_wall")
p.plot_static_wall()
p.animate_1D("rho")
#p.animate("p")
plt.show()
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})