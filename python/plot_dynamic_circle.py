from plotting_utilities import *


#Read header
output_folder = "output_dynamic_circle"


p = Plotter(output_folder)
#p.debug_points(232)
#p.debug_animation()
#p.plot_steady_state("p")
p.animate("p")
plt.show()