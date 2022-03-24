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
#p.debug_points(10)
for n in range(1400,1700+1):
    if n%60 == 0:
        plt.figure()
        p.plot_piston_fsi("p",n,0,10e5)

p.animate_piston_fsi("p",0,10e5)
#p.animate("p",0,10e5)
p.plot_convergence()
plt.show()