

import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
from riemann_exact import riemann_exact
from plotting_utilities import *

p = Plotter("output_riemann_sod")
p.plot_riemann_problem()
p.plot_convergence()
plt.show()
