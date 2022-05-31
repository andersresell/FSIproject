import sys

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from analytical_solutions import *
from decimal import Decimal
import sys
from matplotlib import ticker

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 25,
    #"axes.facecolor": "xkcd:mint green"
    #"axes.facecolor": "lightsteelblue"
    #"lines.linewidth": 2,'
    "axes.labelsize": 30,
    "legend.fontsize": 21,

    "font.size": 20,
    #"axes.facecolor": "xkcd:mint green"
    #"axes.facecolor": "lightsteelblue"
    #"lines.linewidth": 2,'
    "axes.labelsize": 25,
    "legend.fontsize": 21,


    #"figure.figsize": (6,6)

    })

class Plotter:
    def __init__(self, output_folder):
        self.output_folder = output_folder
        #read fvm header:
        data = genfromtxt("output_folders/"+output_folder+"/fvm_header.csv",comments = "#", delimiter=',')
        self.ni = int(data[0])
        self.nj = int(data[1])
        self.L_x = data[2]
        self.L_y = data[3]
        self.write_stride = int(data[4])
        self.n_timesteps = int(data[5])
        self.t_end = data[6]
        self.dx = self.L_x/self.ni
        self.dy = self.L_y/self.nj
        self.x = np.linspace(self.dx/2,self.L_x-self.dx/2,self.ni)
        self.y = np.linspace(self.dy/2,self.L_y-self.dy/2,self.nj)
        #read fsi header:
        data = genfromtxt("output_folders/"+output_folder+"/fsi_header.csv",comments = "#", delimiter=',')
        self.n_static_solids = int(data[0])
        self.n_movable_solids = int(data[1])

    def contour_animate(self,datatype="p",bottom_level=0, top_level=0) :
        plt.figure()
        for n in range(0,self.n_timesteps+1):
            if n % self.write_stride == 0:
                self.contour(datatype=datatype,n=n,bottom_level=bottom_level,top_level=top_level,high_res=False,show_time=True)
                plt.pause(0.00001)
            print(n)
    def contour_plot(self,datatype="p", n=-1, t0 = -1, bottom_level=0, top_level=0,show_time=False):
        plt.figure(figsize=(7,6))
        if n == -1 and t0== -1:
            n = self.n_timesteps
        elif t0 >= 0:
            for n in range(0,self.n_timesteps+1):
                if n % self.write_stride == 0:
                    if (self.timestep2time(n) >= t0):
                        self.contour(datatype=datatype,n=n,bottom_level=bottom_level, top_level=top_level)
                        return
        self.contour(datatype=datatype,n=n,bottom_level=bottom_level, top_level=top_level,high_res=True,show_time=show_time)


    def contour(self,datatype,n,bottom_level=0,top_level=0,high_res=False, show_time=False):
        n_levels = 100
        if high_res:
            n_levels = 300
        if (bottom_level == 0 and top_level == 0):
            auto_level = True
        else:
            auto_level = False
            levels = np.linspace(bottom_level,top_level,n_levels)
        data = self.extract_data(datatype,n)
        plt.clf()
        if datatype=="M":
            cs = plt.contourf(self.x,self.y,data,n_levels,cmap=plt.get_cmap('jet'))
            #plt.contour(self.x,self.y,data,30, linewidths=0.5,colors='k')
        else:
            if auto_level:
                cs = plt.contourf(self.x,self.y,data,n_levels,cmap=plt.get_cmap('hot'))
            else:
                cs = plt.contourf(self.x,self.y,data, levels=levels, cmap=plt.get_cmap('hot'))
        cb = plt.colorbar(cs)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        if not auto_level:
            plt.clim(bottom_level,top_level)
        if datatype == "M":
            cb.set_label("$M$")
        elif datatype == "p":
            cb.set_label(r"""$p\;[\textrm{Pa}]$""")
        #plt.axis('equal')
        ax = plt.gca()
        ax.set_aspect('equal', 'box')
        plt.xlabel(r"""$x\;[\textrm{m}]$""")
        plt.ylabel(r"""$y\;[\textrm{m}]$""")
        if show_time:
            self.time_title(n,"ms")
        ax = plt.gca()
        self.plot_solids(n)
        #ax.set_facecolor("lightsteelblue")
        plt.tight_layout()


    def extract_data(self,datatype,n):
        data = genfromtxt("output_folders/"+self.output_folder+"/fvm_output_t"+str(n)+".csv",comments = "#", delimiter=',')
        if (datatype == "M"):
            M = np.sqrt((data[:,1]**2 + data[:,2]**2)/(1.4*data[:,3]/data[:,0]))
            return np.transpose(M.reshape((self.ni,self.nj)))
        elif (datatype == "rho"):
            rho = data[:,0]
            return np.transpose(rho.reshape((self.ni,self.nj)))
        elif (datatype == "u"):
            u = data[:,1]
            return np.transpose(u.reshape((self.ni,self.nj)))
        elif (datatype == "p"):
            p = data[:,3]
            return np.transpose(p.reshape((self.ni,self.nj)))


    def plot_convergence(self):
        data = genfromtxt("output_folders/"+self.output_folder+"/fvm_convergence_history.csv",comments = "#", delimiter=',')
        #plt.figure(figsize=(7,6))
        plt.figure()
        plt.plot(data[:,0],data[:,1],"k")
        plt.xlabel("$n$")
        plt.ylabel(r'$||\rho^{n+1} -\rho^n||_2$')
        plt.semilogy()
        plt.tight_layout()

        #plt.figure(figsize=(5,5))
        plt.figure()
        plt.plot(data[:,0],data[:,2],"k",linewidth=1)

        plt.xlabel(r"""$n$""")
        plt.ylabel(r"""$m \; [\textrm{kg}]$""")
        plt.tight_layout()
        plt.locator_params(axis="x", nbins=5)
        plt.locator_params(axis="y", nbins=5)
        mass0 = data[0,2]
        mass_end = data[-1,2]
        error = 100*abs(mass_end-mass0)/mass0
        print(mass0,mass_end,error)

    def plot_solids(self, n):
        plt.autoscale(False)
        for i in range(0,self.n_static_solids+self.n_movable_solids):
            if i < self.n_static_solids:
                data = genfromtxt("output_folders/"+self.output_folder+"/static_boundary"+str(i)+".csv",comments = "#", delimiter=',')
            else:
                data = genfromtxt("output_folders/"+self.output_folder+"/movable_boundary"+str(i-self.n_static_solids)+"_t"+str(n)+".csv",comments = "#", delimiter=',')
            plt.fill(data[:,0],data[:,1],"silver")
            plt.plot(np.append(data[:,0],data[0,0]),np.append(data[:,1],data[0,1]),"black")
        #plt.autoscale(True)

    def plot_riemann_problem(self):

        data = genfromtxt("output_folders/"+self.output_folder+"/fvm_output_t"+str(self.n_timesteps)+".csv",comments = "#", delimiter=',')
        rho = data[:,0]
        rho = np.transpose(rho.reshape((self.ni,self.nj)))
        u = data[:,1]
        u = np.transpose(u.reshape((self.ni,self.nj)))
        p = data[:,3]
        p = np.transpose(p.reshape((self.ni,self.nj)))

        rho_l = rho[(int(self.nj/2)),0]
        rho_r = rho[int(self.nj/2),self.ni-1]
        u_l = u[int(self.nj/2),0]
        u_r = u[int(self.nj/2),self.ni-1]
        p_l = p[int(self.nj/2),0]
        p_r = p[int(self.nj/2),self.ni-1]

        rho_exact= np.zeros(self.x.shape)
        u_exact=rho_exact
        p_exact=rho_exact
        rho_exact,u_exact,p_exact = riemann_exact_solution(rho_l,u_l,p_l,rho_r,u_r,p_r,self.ni,self.L_x,self.t_end)
        #plt.figure()
        #s_exact = np.log(p_exact/rho_exact**1.4)
        #plt.plot(self.x,s_exact)


        plt.figure()
        plt.plot(self.x,rho[int(self.nj/2),:],'k.')
        #plt.plot(self.x,rho_exact,'--r')
        plt.legend(['Numerical','Exact'])
        plt.xlabel("$x$")
        plt.ylabel(r'$\rho$')

        plt.figure()
        plt.plot(self.x,u[int(self.nj/2),:],'k.')
        plt.plot(self.x,u_exact,'--r')
        plt.legend(['Numerical','Exact'])
        plt.xlabel("$x$")
        plt.ylabel("$u$")

        plt.figure()
        plt.plot(self.x,p[int(self.nj/2),:],'k.')
        plt.plot(self.x,p_exact,'--r')
        plt.legend(['Numerical','Exact'])
        plt.xlabel("$x$")
        plt.ylabel("$p$")


        plt.rcParams.update({
            "font.size": 20,
            #"axes.facecolor": "xkcd:mint green"
            #"axes.facecolor": "lightsteelblue"
            #"lines.linewidth": 2,'
            #"axes.labelsize": 28,
            #"legend.fontsize": 20

        })
        fig = plt.figure(figsize=(5,9))
        ax1 = plt.subplot(3,1,1)
        a = plt.plot(self.x,rho[int(self.nj/2),:],'k.',label='Numerical')
        b = plt.plot(self.x,rho_exact,'--r',label='Exact')
        #plt.legend(['Numerical','Exact'])
        plt.xlabel("$x$")
        plt.ylabel(r'$\rho$')

        lines, labels = fig.axes[-1].get_legend_handles_labels()
        #fig.legend(lines, labels, loc = 'upper right',)
        fig.legend(lines,labels,bbox_to_anchor=(1,0.887), loc="lower right",  bbox_transform=fig.transFigure)

        ax2 = plt.subplot(3,1,2,sharex=ax1)
        plt.plot(self.x,u[int(self.nj/2),:],'k.')
        plt.plot(self.x,u_exact,'--r')
        plt.xlabel("$x$")
        plt.ylabel("$u$")
        #ax1.set_xticks([1,2.5,2])

        ax3 = plt.subplot(3,1,3,sharex=ax2)
        plt.plot(self.x,p[int(self.nj/2),:],'k.')
        plt.plot(self.x,p_exact,'--r')
        plt.xlabel("$x$")
        plt.ylabel("$p$")

        plt.suptitle("$NI = "+str(self.ni)+"$",y=0.96,x=0.25)
        plt.tight_layout()

    def plot_1D(self,datatype,bottom_level=0,top_level=0,n=-1,t_goal=0,show_timestep=True):
        if n == -1 and t_goal == 0:
            n = self.n_timesteps
        elif t_goal>0:
            n = self.time2timestep(t_goal)
        print(n,t_goal)
        data = self.extract_data(datatype,n)
        plt.plot(self.x,data[(int(self.nj/2)),:],'k')
        plt.xlabel("x [m]")
        if datatype == "p":
            plt.ylabel("p [Pa]")
        if top_level != 0:
            plt.ylim(bottom_level,top_level)
        if show_timestep:
            self.time_title(n,"ms")
        plt.tight_layout()
        plt.locator_params(axis="x", nbins=5)
        plt.locator_params(axis="y", nbins=4)


    def animate_1D(self,datatype,bottom_level=0,top_level=0,show_timestep=True):
        plt.figure()
        for n in range(0,self.n_timesteps+1):
            if n % self.write_stride == 0:
                plt.cla()
                self.plot_1D(datatype=datatype,bottom_level=bottom_level,top_level=top_level,n=n,show_timestep=show_timestep)
                plt.pause(0.1)


    def plot_static_wall(self):
        rho_l = self.extract_data("rho",0)
        rho_l = rho_l[0,0]
        u_l = self.extract_data("u",0)
        u_l = u_l[0,0]
        p_l = self.extract_data("p",0)
        p_l = p_l[0,0]

        rho = self.extract_data("rho",self.n_timesteps)
        u = self.extract_data("u",self.n_timesteps)
        p = self.extract_data("p",self.n_timesteps)

        rho_l_star, p_star, s = static_wall_exact(rho_l,u_l,p_l)
        x_s = self.L_x + s*self.t_end
        ind = np.where(self.x >= x_s)
        rho_exact = np.ones(self.x.shape)*rho_l
        u_exact = np.ones(self.x.shape)*u_l
        p_exact = np.ones(self.x.shape)*p_l
        rho_exact[ind] = rho_l_star
        u_exact[ind] = 0
        p_exact[ind] = p_star

        plt.figure()
        plt.plot(self.x,rho[int(self.nj/2),:],'k.')
        plt.plot(self.x,rho_exact,'--r')
        plt.legend(['numerical','exact'])
        plt.xlabel("x")
        plt.ylabel("Density")

        plt.figure()
        plt.plot(self.x,u[int(self.nj/2),:],'k.')
        plt.plot(self.x,u_exact,'--r')
        plt.legend(['numerical','exact'])
        plt.xlabel("x")
        plt.ylabel("Velocity")

        plt.figure()
        plt.plot(self.x,p[int(self.nj/2),:],'k.')
        plt.plot(self.x,p_exact,'--r')
        plt.legend(['numerical','exact'])
        plt.xlabel("x")
        plt.ylabel("Pressure")

    def plot_piston_prescribed(self,vel):
        plt.rcParams.update({
            "font.size": 25,
            "axes.labelsize": 30,
            "legend.fontsize": 21,
        })
        rho_c = self.extract_data("rho",0)
        rho_c = rho_c[0,0]
        u_c = self.extract_data("u",0)
        u_c = u_c[0,0]
        p_c = self.extract_data("p",0)
        p_c = p_c[0,0]

        #a = 0.5*(self.L_x-width) + vel*self.t_end
        #b = 0.5*(self.L_x+width) + vel*self.t_end
        piston = genfromtxt("output_folders/"+self.output_folder+"/movable_boundary0_t"+str(self.n_timesteps)+".csv",comments = "#", delimiter=',')
        a = piston[0,0]
        b = piston[1,0]
        width = b-a

        def get_data(datatype):
            ind = np.where(np.logical_and(self.x > a, self.x < b))
            data = self.extract_data(datatype,self.n_timesteps)
            data = data[int(self.nj/2),:]
            data[ind] = float("nan")
            return data

        def visualize_piston(data):
            MAX = max(data)
            MIN = min(data)
            delta = (MAX-MIN)*0.1
            MAX+=delta
            MIN-=delta
            x = np.array([a,a])
            y = np.array([MIN,MAX])
            #plt.plot(x,y,'grey',linewidth=3)
            x = np.array([b,b])
            y = np.array([MIN,MAX])
            #plt.plot(x,y,'grey',linewidth=3)

            x = np.array([a,b,b,a])
            y = np.array([MIN,MIN,MAX,MAX])
            plt.fill(x,y,"silver")
            plt.plot(np.append(x,x[0]),np.append(y,y[0]),"black")
            plt.ylim(MIN,MAX)

        rho_exact,u_exact,p_exact = piston_exact(rho_c,p_c,rho_c,p_c,vel,self.ni,self.L_x,width,self.t_end)

        rho = get_data("rho")
        plt.figure()
        plt.plot(self.x,rho,'k.')
        plt.plot(self.x,rho_exact,'--r',linewidth=2)
        plt.legend(['Numerical','Exact'])
        visualize_piston(rho)
        self.set_unit_labels_SI("rho")
        plt.tight_layout()
        ax = plt.gca()
        ax.set_yticks([0,1.5,3])

        u = get_data("u")
        plt.figure()
        plt.plot(self.x,u,'k.')
        plt.plot(self.x,u_exact,'--r',linewidth=2)
        plt.legend(['Numerical','Exact'])
        visualize_piston(u)
        self.set_unit_labels_SI("u")
        plt.ylabel(r"""$u\;[\textrm{m/s}]$""")
        plt.tight_layout()
        ax = plt.gca()
        ax.set_yticks([0,250,500])

        p = get_data("p")
        plt.figure()
        plt.plot(self.x,p,'k.')
        plt.plot(self.x,p_exact,'--r',linewidth=2)
        plt.legend(['Numerical','Exact'])
        visualize_piston(p)
        self.set_unit_labels_SI("p")
        ax = plt.gca()
        ax.set_yticks([0,250000,500000])
        plt.tight_layout()

    def animate_piston_fsi(self,datatype,bottom_level=0, top_level=0):
        plt.figure()
        for n in range(0,self.n_timesteps+1):
            if n % self.write_stride == 0and n%100==0:
                print(n)
                plt.clf()
                self.piston_fsi(datatype,n,bottom_level,top_level)
                plt.pause(0.1)
    def plot_piston_fsi(self,datatype,t0=-1,bottom_level=0,top_level=0):
        plt.figure()
        if t0 >= 0:
            for n in range(0,self.n_timesteps+1):
                if n % self.write_stride == 0:
                    if self.timestep2time(n) >= t0:
                        self.piston_fsi(datatype=datatype,n=n,bottom_level=bottom_level,top_level=top_level)
                        return
        self.piston_fsi(datatype,self.n_timesteps,bottom_level,top_level)
    def piston_fsi_extract_data(self,t_a=0,t_b=0,get_c_wall=False,stride=1):
        if t_b == 0:
            last_step = self.n_timesteps
        else:
            last_step = self.time2timestep(t_b)
        t = np.array([])
        a = np.array([])
        p_wall = np.array([])
        u_wall = np.array([])
        c_wall = np.array([])
        for n in range(0,self.n_timesteps):
            if n%self.write_stride == 0 and n%stride==0:
                if self.timestep2time(n) >= t_a and n <= last_step:
                    t = np.append(t,self.timestep2time(n))
                    #t[n] = self.timestep2time(n)
                    print("t = "+str(t[-1])+" step n = "+str(n)+" of "+str(self.n_timesteps))
                    if self.n_movable_solids == 1:
                        piston = genfromtxt("output_folders/"+self.output_folder+"/movable_boundary0_t"+str(n)+".csv",comments = "#", delimiter=',')
                        piston_vel = genfromtxt("output_folders/"+self.output_folder+"/movable_solid_body_CM_velocity0_t"+str(n)+".csv",comments = "#", delimiter=',')
                    elif self.n_static_solids == 1:
                        piston = genfromtxt("output_folders/"+self.output_folder+"/static_boundary0.csv",comments = "#", delimiter=',')
                        piston_vel = np.array([0])
                    else:
                        sys.exit("error! There has to be exactly on solid object for this simulation")
                    #a[n] = piston[0,0]
                    a = np.append(a,piston[0,0])
                    #u_wall[n] = piston_vel[0]

                    u_wall = np.append(u_wall,piston_vel[0])
                    b = piston[1,0]
                    #if t[n] >= 0.02:# and t[n] <= 0.03:
                    #p_wall[n] = self.probe_1D(datatype="p",x=a[n],n=n)
                    p_wall = np.append(p_wall,self.probe_1D(datatype="p",x=a[-1],n=n))
                    if get_c_wall:
                        rho_wall = self.probe_1D(datatype="p",x=a[-1],n=n)
                        c_wall = np.append(c_wall,np.sqrt(1.4*p_wall[-1]/rho_wall))
                    #u_wall[n] = self.probe_1D(datatype="u",x=a[n],n=n)
        return t,a,p_wall,u_wall,c_wall

    def piston_fsi_comparison(self,updated_comparison=False,stride=1):
        t,a,p_wall,u_wall,c_wall = self.piston_fsi_extract_data(stride=stride)
        ind = np.where(u_wall>1)
        g = 1.4
        if not updated_comparison:
            p_r = 701258.5
            rho_r = 4.572815
            c_r = np.sqrt(g*p_r/rho_r)
            pressure_ratio_analytical = (1 - (g-1)/2*u_wall/c_r)**((2*g)/(g-1))
            pressure_ratio_sim = p_wall/p_r
        else:
            p_stat = Plotter("piston_fsi_static")
            t_stat,a_stat,p_wall_stat,u_wall_stat,c_wall_stat = p_stat.piston_fsi_extract_data()
            p_r = np.zeros(t.shape)
            c_r = np.zeros(t.shape)
            #for i in range(0,t.size):
            #    p_r[i] = np.argmax()
            #pressure_ratio_analytical = (1 - (g-1)/2*u_wall/c_r)**((2*g)/(g-1))
            pressure_ratio_sim = p_wall/p_r


        plt.figure()
        plt.plot(t*1000,u_wall,'k')
        plt.xlabel(r"""$t\;[\textrm{ms}]$""")
        plt.ylabel(r"""$u_{wall}\;[\textrm{m/s}]$""")
        plt.xlim(20,40)
        plt.tight_layout()

        plt.figure()
        plt.plot(t*1000,pressure_ratio_sim,'k')
        plt.plot(t[ind]*1000,pressure_ratio_analytical[ind],'--r')
        plt.xlim(20,40)
        plt.xlabel(r"""$t\;[\textrm{ms}]$""")
        #plt.ylabel(r"""$p'_r/p_r$""")
        #plt.legend(['Simulation','Analytical Formula'])
        plt.legend([r"""$p'_{r,sim}/p_r$""",r"""$\left(1-\frac{\gamma-1}{2}\frac{u_{wall}}{c_{r}}\right)^{\frac{2\gamma}{\gamma-1}}$"""])
        plt.tight_layout()

    def piston_fsi(self,datatype,n,bottom_level=0,top_level=0):
        autolevel = False
        if bottom_level==0 and top_level==0:
            autolevel = True
        if self.n_movable_solids == 1:
            piston = genfromtxt("output_folders/"+self.output_folder+"/movable_boundary0_t"+str(n)+".csv",comments = "#", delimiter=',')
        elif self.n_static_solids == 1:
            piston = genfromtxt("output_folders/"+self.output_folder+"/static_boundary0.csv",comments = "#", delimiter=',')
        else:
            sys.exit("error! There has to be exactly on solid object for this simulation")
        a = piston[0,0]
        b = piston[1,0]
        def get_data(datatype):
            ind = np.where(np.logical_and(self.x > a, self.x < b))
            data = self.extract_data(datatype,n)
            data = data[int(self.nj/2),:]
            data[ind] = float("nan")
            return data

        def visualize_piston(data):
            if bottom_level==0 and top_level == 0:
                MAX = max(data)
                MIN = min(data)
            else:
                MAX = top_level
                MIN = bottom_level
            delta = (MAX-MIN)*0.1
            MAX+=delta
            MIN-=delta
            x = np.array([a,b,b,a])
            y = np.array([MIN,MIN,MAX,MAX])
            plt.fill(x,y,"silver")
            plt.plot(np.append(x,x[0]),np.append(y,y[0]),"black")
            plt.ylim(MIN,MAX)

        data = get_data(datatype)
        plt.plot(self.x,data,'k.',markersize=2.5)
        visualize_piston(data)

        #plt.ticklabel_format(scilimits=(-3,9))
        #plt.xlabel(r"""$x\;[\textrm{m}]$""")
        #plt.ylabel(r"""$p\;[\textrm{Pa}]$""")
        self.set_unit_labels_SI("p")
        self.time_title(n,"ms")
        plt.tight_layout()

    def debug_points_animation(self):
        plt.figure()
        for n in range(0,self.n_timesteps):
            if n % self.write_stride == 0:
                plt.clf()
                self.debug_points(n)
                plt.pause(0.1)

    def debug_points_last(self):
        plt.figure()
        self.debug_points(self.n_timesteps)

    def debug_points(self, n):
        #plotting points
        data = genfromtxt("output_folders/"+self.output_folder+"/debug_nodes_t"+str(n)+".csv",comments = "#", delimiter=',')
        cell_type = data[:,0]
        x = data[:,1]
        y = data[:,2]
        for i in range(0,cell_type.shape[0]):
            if cell_type[i] == 0:
                #fluid
                plt.plot(x[i],y[i],'ob')
            elif cell_type[i] == 1:
                #ghost
                plt.plot(x[i],y[i],'om')
            elif cell_type[i] == 2:
                #solid
                plt.plot(x[i],y[i],'or')
        #plotting boundaries
        for i in range(0,self.n_static_solids+self.n_movable_solids):
            if i < self.n_static_solids:
                data = genfromtxt("output_folders/"+self.output_folder+"/static_boundary"+str(i)+".csv",comments = "#", delimiter=',')
            else:
                data = genfromtxt("output_folders/"+self.output_folder+"/movable_boundary"+str(i-self.n_static_solids)+"_t"+str(n)+".csv",comments = "#", delimiter=',')
            x_b = data[:,0]
            y_b = data[:,1]
            plt.plot(x_b,y_b,'.-k')
            plt.plot(np.array([x_b[-1],x_b[0]]),np.array([y_b[-1],y_b[0]]),'.-k')
        #plotting intercepts
        data = genfromtxt("output_folders/"+self.output_folder+"/debug_intercepts_t"+str(n)+".csv",comments = "#", delimiter=',')
        if not data.size == 0:
            x_i = data[:,0]
            y_i = data[:,1]
            plt.plot(x_i,y_i,'y*')
            plt.axis('equal')
        #plotting fresh points
        #if n != 0:
         #   data = genfromtxt("output_folders/"+self.output_folder+"/debug_fresh_points_t"+str(n)+".csv",comments = "#", delimiter=',')

          #  if data.size == 1:
           #     x_i = data[0]
            #    y_i = data[1]
            #elif data.size(0) > 1:
            #    x_i = data[:,0]
            #    y_i = data[:,1]
            #if data.size(0) > 0:
            #    plt.plot(x_i,y_i,'kx')
            #    plt.axis('equal')
    #def interpolate(self,x,y,datatype,n):
     #   data = self.extract_data(datatype,n)
      #  im = int(x/self.dx-0.5)
       # jm = int(y/self.dy-0.5)


    def probe_1D(self,datatype,x,n=-1):
        #the method is meant for 1D interpolation, and should be modified before being used for 2D interpolation
        if n==-1:
            n = self.n_timesteps
        data = self.extract_data(datatype,n)
        i = int(x/self.dx-0.5)
        a = self.dx*(i+0.5)
        b = a + self.dx
        epsilon = (x-a)/(b-a)
        j = 0
        return data[j,i]*(1-epsilon) + data[j,i+1]*epsilon
    def time_history_1D(self,x_point,datatype="p",t_a=0,t_b=0,stride=1):
        if t_b != 0:
            n_first = self.time2timestep(t_a)
            n_last = self.time2timestep(t_b)
        else:
            n_first = 0
            n_last = self.n_timesteps
        val = np.array([])
        t = np.array([])
        for n in range(n_first,n_last):
                if n%self.write_stride == 0 and n%stride==0:
                    print("probing time history "+str(n)+" of "+str(self.n_timesteps))
                    val = np.append(val,self.probe_1D(datatype,x_point,n))
                    t = np.append(t,self.timestep2time(n))
        return t,val
    def timestep2time(self,n):
        f = open("output_folders/"+self.output_folder+"/fvm_output_t"+str(n)+".csv","r")
        first_line = f.readline()
        return float(first_line[3:])
    def time2timestep(self,t_goal):
        for n in range(0,self.n_timesteps):
            if n%self.write_stride == 0:
                if self.timestep2time(n) >= t_goal:
                    return n
        sys.exit("error! t_goal larger than max t. no n found")
    def time_title(self,n,unit="s"):
        if unit == "s":
            plt.title("t = "+f"{self.timestep2time(n):.2f}"+" [s]")
        elif unit == "ms":
            #plt.title("time = "+f"{self.timestep2time(n)*1e3:.2f}"+" [ms]")
            plt.title(r"""$$ t = """+f"{self.timestep2time(n)*1e3:.2f}"+r""" \;\textrm{ms}$$""")
    def set_unit_labels_SI(self,datatype):
        plt.xlabel(r"""$x\;[\textrm{m}]$""")
        if datatype == "rho":
            plt.ylabel(r"""$\rho\;[\textrm{kg/m}^3]$""")
        elif datatype == "u":
            plt.ylabel(r"""$u\;[\textrm{m/s}]$""")
        elif datatype == "p":
            plt.ylabel(r"""$p\;[\textrm{Pa}]$""")
    def set_ticks(self, n_x,n_y):
        ax = plt.gca()
        ax.set_xticks(np.linspace(0,self.L_x,n_x))
        ax.set_yticks(np.linspace(0,self.L_y,n_x))
    def schlieren_last(self):
        plt.figure()
        self.schlieren(self.n_timesteps)
    def schlieren_animate(self):
        plt.rcParams.update({
            "axes.facecolor": "black"
        })
        fig = plt.figure()
        fig.patch.set_facecolor("black")
        for n in range(0,self.n_timesteps):
            if n % self.write_stride == 0:
                plt.cla()
                self.schlieren(n,show_timestep=True)
                plt.pause(0.00001)
    def schlieren_plot_all(self):
        for n in range(800,self.n_timesteps):
            if n % self.write_stride == 0:
                plt.figure()
                self.schlieren(n)
        plt.show()
    def schlieren_plot_time_interval(self,t0,show_timestep):
        t = 0
        delta_t = 1/37000
        first_found = False
        #print('%.2E' % Decimal('40800000000.00000000000000'))
        for n in range(0,self.n_timesteps):
            if n % self.write_stride == 0:
                #print(self.timestep2time(n))
                if self.timestep2time(n) >= t0+t:
                    plt.figure()
                    plt.rcParams.update({
                        "axes.facecolor": "black"
                    })
                    if not first_found:
                        t0 = self.timestep2time(n)
                        first_found = True
                    self.schlieren(n,show_timestep)
                    plt.title("t = "+f"{self.timestep2time(n)-t0:.6f}")
                    t += delta_t
        plt.show()

    def schlieren(self,n,show_timestep=False):
        data = self.extract_data("rho",n)
        x_center = 0.26/2
        y_center = 0.26/2
        w = 0.120
        h = w
        a = x_center-w/2
        b = x_center+w/2
        ii = np.where(np.logical_and(self.x>=x_center-w/2,self.x<=x_center+w/2))
        jj = np.where(np.logical_and(self.y>=y_center-h/2,self.y<=y_center+h/2))
        data = data[np.ix_(jj[0],ii[0])]
        d_rho_dx = np.zeros(data.shape)
        d_rho_dy = np.zeros(data.shape)
        d_rho_dx[1:-1,:] = (data[2::,:] - data[0:-2,:])/(2*self.dx)
        d_rho_dx[0,:] = (data[1,:] - data[0,:])/self.dx
        d_rho_dx[-1,:] = (data[-1,:] - data[-2,:])/self.dx

        d_rho_dy[:,1:-1] = (data[:,2::]- data[:,0:-2])/(2*self.dy)
        d_rho_dy[:,0] = (data[:,1] - data[:,0])/self.dy
        d_rho_dy[:,-1] = (data[:,-1] - data[:,-2])/self.dy

        d_rho_norm = np.sqrt(d_rho_dx**2+d_rho_dy**2)
        levels = np.linspace(0,200,50)

        cs = plt.contourf(self.x[ii],self.y[jj],d_rho_norm,levels=levels,cmap=plt.get_cmap('binary'))
        #cs = plt.contourf(self.x[ii],self.y[jj],d_rho_norm,100,cmap=plt.get_cmap('binary'))
        self.plot_solids(self.n_timesteps)

        #plt.axis('equal')
        ax = plt.gca()
        ax.set_aspect('equal', 'box')
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        plt.tick_params(left = False)
        plt.tick_params(bottom = False)
        plt.tight_layout()
        plt.rcParams.update({
            "axes.facecolor": "black"
        })
        if show_timestep:
            plt.title("t = "+str(self.timestep2time(n)))
