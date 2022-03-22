from plotting_utilities import *


#Read header
output_folder = "output_tube"



p = Plotter(output_folder)
#p.debug_points_last()
#p.debug_animation()
p.plot_steady_state("M")
#p.animate("p",0.5e5,10e5)
plt.show()