from plotting_utilities import *


#Read header
output_folder = "output_dynamic_unconstrained"



p = Plotter(output_folder)
#p.debug_points_last()
#p.debug_animation()
#p.plot_steady_state("p")
p.animate("p",2e4,7e5)
p.plot_convergence()
plt.show()