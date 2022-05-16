from plotting_utilities import *


#Read header
output_folder = "output_dynamic_unconstrained"



p = Plotter(output_folder)
#p.debug_points_last()
#p.debug_animation()
#p.plot_steady_state("p")
times = np.linspace(0,0.02,6)
times = np.array([2,5,12,25])*1e-3
for i in range(0,times.size):
    p.contour_plot(datatype="p",t0=times[i],bottom_level=0,top_level=5e5)

#p.contour_animate("p",1e4,5e5)
#p.contour_animate("p")
#p.plot_convergence()
#plt.figure()
#p.contour("p",2000,bottom_level=0,top_level=900000)
#p.contour_plot("p",bottom_level=0,top_level=10e5)
plt.show()