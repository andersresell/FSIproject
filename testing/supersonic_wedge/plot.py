# import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import genfromtxt
import os
import sys

root_dir = os.getenv("FSI_PROJECT_HOME")
python_dir = os.path.join(root_dir, "python")
print("dir", python_dir)
sys.path.append(python_dir)
from analytical_solutions import riemann_exact
from plotting_utilities import *

# need this to make ctrl+c close plot window
import signal

signal.signal(signal.SIGINT, signal.SIG_DFL)


# Read header
sim_dir = os.getcwd()
# output_folder = "wedge_verification_detached"

p = Plotter(sim_dir)
p.contour_animate("p")
# p.contour_animate("p")
p.contour_plot(datatype="M")
p.plot_convergence()
# print(p.probe("M",9,14))
# print(p.probe("M",4,7))
# print(p.probe("M",4.67,6.23))
plt.show()
