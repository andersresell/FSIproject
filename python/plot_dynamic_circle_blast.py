from plotting_utilities import *


#Read header
output_folder = "output_dynamic_circle_blast"


p = Plotter(output_folder)
#p.debug_points(1)
#plt.figure()
#p.debug_points(p.n_timesteps)
#plt.figure()
#p.debug_points(n)
#p.debug_animation()
#p.plot_steady_state("p")
p.contour_animate("p",0,3)
p.plot_convergence()
plt.show()