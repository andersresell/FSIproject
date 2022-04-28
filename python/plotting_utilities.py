import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from analytical_solutions import *

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

    def contour_animate(self,datatype,bottom_level=0, top_level=0) :
        plt.figure()
        for n in range(0,self.n_timesteps+1):
            if n % self.write_stride == 0:
                self.contour( datatype,n,bottom_level,top_level)
                plt.pause(0.00001)
            print(n)
    def contour_plot(self,datatype, n=-1,contour_lines=False):
        plt.figure(figsize=(7,6))
        if n == -1:
            n = self.n_timesteps
        self.contour(datatype,n,0,0,True,contour_lines)

    def contour(self,datatype,n,bottom_level=0,top_level=0,high_res=False,contour_lines=False):
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
        if contour_lines:
            cs = plt.contourf(self.x,self.y,data,300,cmap=plt.get_cmap('jet'))
            #plt.contour(self.x,self.y,data,30, linewidths=0.5,colors='k')
        else:
            if auto_level:
                cs = plt.contourf(self.x,self.y,data,n_levels,cmap=plt.get_cmap('hot'))
            else:
                cs = plt.contourf(self.x,self.y,data, levels=levels, cmap=plt.get_cmap('hot'))
        cb = plt.colorbar(cs)
        if datatype == "M":
            cb.set_label("$M$")
        elif datatype == "p":
            cb.set_label("$p\;[Pa]$")
        self.plot_solids(n)
        plt.axis('equal')
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.title("t = "+str(self.timestep2time(n))+"[s]")
        ax = plt.gca()
        #ax.set_facecolor("lightsteelblue")
        #plt.tight_layout()


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
        plt.xlabel("$n$")
        plt.ylabel("$m \; [kg]$")
        plt.tight_layout()

    def plot_solids(self, n):
        for i in range(0,self.n_static_solids+self.n_movable_solids):
            if i < self.n_static_solids:
                data = genfromtxt("output_folders/"+self.output_folder+"/static_boundary"+str(i)+".csv",comments = "#", delimiter=',')
            else:
                data = genfromtxt("output_folders/"+self.output_folder+"/movable_boundary"+str(i-self.n_static_solids)+"_t"+str(n)+".csv",comments = "#", delimiter=',')
            plt.fill(data[:,0],data[:,1],"silver")
            plt.plot(np.append(data[:,0],data[0,0]),np.append(data[:,1],data[0,1]),"black")

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

        rho_exact,u_exact,p_exact = riemann_exact_solution(rho_l,u_l,p_l,rho_r,u_r,p_r,self.ni,self.L_x,self.t_end)

        #plt.figure()
        #s_exact = np.log(p_exact/rho_exact**1.4)
        #plt.plot(self.x,s_exact)


        plt.figure()
        plt.plot(self.x,rho[int(self.nj/2),:],'k.')
        plt.plot(self.x,rho_exact,'--r')
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

        ax3 = plt.subplot(3,1,3,sharex=ax2)
        plt.plot(self.x,p[int(self.nj/2),:],'k.')
        plt.plot(self.x,p_exact,'--r')
        plt.xlabel("$x$")
        plt.ylabel("$p$")

        plt.suptitle("$ni = "+str(self.ni)+"$",y=0.96,x=0.25)
        plt.tight_layout()


    def animate_1D(self,datatype,bottom_level=0,top_level=0):
        plt.figure()
        for n in range(0,self.n_timesteps+1):
            if n % self.write_stride == 0:
                data = self.extract_data(datatype,n)
                plt.clf()
                plt.plot(self.x,data[(int(self.nj/2)),:])
                plt.xlabel("x")
                plt.ylabel(datatype)
                print(top_level)
                if top_level != 0:
                    plt.ylim(bottom_level,top_level)
                plt.tight_layout()
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
        plt.xlabel("$x$")
        plt.ylabel(r'$\rho$')
        #plt.tight_layout()

        u = get_data("u")
        plt.figure()
        plt.plot(self.x,u,'k.')
        plt.plot(self.x,u_exact,'--r',linewidth=2)
        plt.legend(['Numerical','Exact'])
        visualize_piston(u)
        plt.xlabel("$x$")
        plt.ylabel("$u$")
        #plt.tight_layout()

        p = get_data("p")
        plt.figure()
        plt.plot(self.x,p,'k.')
        plt.plot(self.x,p_exact,'--r',linewidth=2)
        plt.legend(['Numerical','Exact'])
        visualize_piston(p)
        plt.xlabel("$x$")
        plt.ylabel("$p$")
        #plt.tight_layout()

    def animate_piston_fsi(self,datatype,bottom_level=0, top_level=0):
        plt.figure()
        for n in range(0,self.n_timesteps+1):
            if n % self.write_stride == 0:
                print(n)
                plt.clf()
                self.plot_piston_fsi(datatype,n,bottom_level,top_level)
                plt.pause(0.1)

    def plot_piston_fsi(self,datatype,n,bottom_level=0,top_level=0):
        autolevel = False
        if bottom_level==0 and top_level==0:
            autolevel = True
        piston = genfromtxt("output_folders/"+self.output_folder+"/movable_boundary0_t"+str(n)+".csv",comments = "#", delimiter=',')
        a = piston[0,0]
        b = piston[1,0]
        def get_data(datatype):
            ind = np.where(np.logical_and(self.x > a, self.x < b))
            data = self.extract_data(datatype,n)
            data = data[int(self.nj/2),:]
            data[ind] = float("nan")
            return data

        def visualize_piston(data):
            MAX = max(data)
            MIN = min(data)
            MAX = 1e6
            MIN = 0
            delta = (MAX-MIN)*0.1
            MAX+=delta
            MIN-=delta
            x = np.array([a,b,b,a])
            y = np.array([MIN,MIN,MAX,MAX])
            plt.fill(x,y,"silver")
            plt.plot(np.append(x,x[0]),np.append(y,y[0]),"black")
            plt.ylim(MIN,MAX)

        data = get_data(datatype)
        plt.plot(self.x,data,'k.',markersize=3.5)
        visualize_piston(data)
        plt.xlabel("x")
        plt.ylabel(datatype)

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


    def probe(self,datatype,x_point,y_point,n=-1):
        if n==-1:
            n = self.n_timesteps
        data = self.extract_data(datatype,n)
        i = int(x_point/self.dx-0.5)
        j = int(y_point/self.dy-0.5)
        return data[j,i]
    def time_history(self,x_point,y_point,datatype="p"):
        val = np.array([])
        t = np.array([])
        for n in range(0,self.n_timesteps):
            if n%self.write_stride == 0:
                val = np.append(val,self.probe(datatype,x_point,y_point,n))
                t = np.append(t,self.timestep2time(n))
        return t,val
    def timestep2time(self,n):
        f = open("output_folders/"+self.output_folder+"/fvm_output_t"+str(n)+".csv","r")
        first_line = f.readline()
        return float(first_line[3:])

    def schlieren_experiment(self):
        plt.figure()
        data = self.extract_data("rho",self.n_timesteps)
        x_center = 5.075
        y_center = 0.29/2
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
        cs = plt.contourf(self.x[ii],self.y[jj],d_rho_norm,300,cmap=plt.get_cmap('binary'))
        self.plot_solids(self.n_timesteps)
        plt.axis('equal')
