import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from riemann_exact import *

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

    def animate(self,datatype):
        for n in range(0,self.n_timesteps+1):
            if n % self.write_stride == 0:
                data = self.extract_data(datatype,n)
                plt.clf()
                cs = plt.contourf(self.x,self.y,data, 100, cmap=plt.get_cmap('hot'))
                cb = plt.colorbar(cs)
                if datatype == "M":
                    cb.set_label("Mach number")
                elif datatype == "p":
                    cb.set_label("Pressure [Pa]")
                self.plot_solids(n)
                plt.axis('equal')
                plt.xlabel('x')
                plt.ylabel('y')
                plt.pause(0.00001)
            print(n)
        plt.show()

    def plot_steady_state(self,datatype):
        data = self.extract_data(datatype,self.n_timesteps)
        self.plot_convergence()
        plt.figure()
        #cs = plt.contour(x,y,M, 70, cmap=plt.get_cmap('jet'),linewidths=0.5)
        cs = plt.contourf(self.x,self.y,data,100,cmap=plt.get_cmap('jet'))
        plt.contour(self.x,self.y,data,30, linewidths=0.5,colors='black')
        cb = plt.colorbar(cs)
        if datatype == "M":
            cb.set_label("Mach number")
        elif datatype == "p":
            cb.set_label("Pressure [Pa]")
        self.plot_solids(self.n_timesteps)
        plt.axis('equal')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def extract_data(self,datatype,n):
        data = genfromtxt("output_folders/"+self.output_folder+"/fvm_output_t"+str(n)+".csv",comments = "#", delimiter=',')
        if (datatype == "M"):
            M = np.sqrt((data[:,1]**2 + data[:,2]**2)/(1.4*data[:,3]/data[:,0]))
            return np.transpose(M.reshape((self.ni,self.nj)))
        elif (datatype == "p"):
            p = data[:,3]
            return np.transpose(p.reshape((self.ni,self.nj)))


    def plot_convergence(self):
        data = genfromtxt("output_folders/"+self.output_folder+"/fvm_convergence_history.csv",comments = "#", delimiter=',')
        plt.figure()
        plt.plot(data[:,0],data[:,1])
        plt.xlabel("n")
        plt.ylabel("L2 norm of density residual")
        plt.semilogy()

    def plot_solids(self, n):
        for i in range(0,self.n_static_solids+self.n_movable_solids):
            if i < self.n_static_solids:
                print("output_folders/"+self.output_folder+"/static_boundary"+str(self.n_static_solids)+".csv")
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

        rho_exact,u_exact,p_exact = riemann_exact(rho_l,u_l,p_l,rho_r,u_r,p_r,self.ni,self.L_x,self.t_end)

        plt.figure()
        plt.plot(self.x,rho[int(self.nj/2),:],'.')
        plt.plot(self.x,rho_exact,"black")
        plt.xlabel("x")
        plt.ylabel("Density")

        plt.figure()
        plt.plot(self.x,u[int(self.nj/2),:],'.')
        plt.plot(self.x,u_exact,"black")
        plt.xlabel("x")
        plt.ylabel("Velocity")

        plt.figure()

        plt.plot(self.x,p[int(self.nj/2),:],'k.')
        plt.plot(self.x,p_exact,'--r')
        plt.legend(['numerical','exact'])

        plt.xlabel("x")
        plt.ylabel("Pressure")

