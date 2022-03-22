
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import optimize as op

#This code solves an exact riemann problem of the euler equations. Subsonic conditions are assumed.
# It implements the procedure described in chapter 14.11 in "Finite Volume Methods for Hyperbolic Problems"

#Sod's problem
#rho_l = 3
#u_l = 0
#p_l = 3
#rho_r = 1
#u_r = 0
#p_r = 1
def riemann_exact(rho_l,u_l,p_l,rho_r,u_r,p_r):
    gamma = 1.4
    beta = (gamma+1)/(gamma-1)
    c_l = math.sqrt(gamma*p_l/rho_l)
    c_r = math.sqrt(gamma*p_r/rho_r)
    def phi_l(p):
        if p <= p_l:
            return u_l + 2*c_l/(gamma-1)*(1 - (p/p_l)**((gamma-1)/(2*gamma)))
        elif p > p_l:
            return u_l + 2*c_l/math.sqrt(2*gamma*(gamma-1))*((1-p/p_l)/math.sqrt(1+beta*p/p_l))

    def phi_r(p):
        if p <= p_r:
            return u_r - 2*c_l/(gamma-1)*(1 - (p/p_r)**((gamma-1)/(2*gamma)))
        elif p > p_r:
            return u_r - 2*c_l/math.sqrt(2*gamma*(gamma-1))*((1-p/p_r)/math.sqrt(1+beta*p/p_r))

    def phi(p):
        return phi_l(p) - phi_r(p)

    p_star0 = 1
    p_star = op.fsolve(phi,p_star0) #Finding the pressure at the intermediate state
    p_star = p_star[0] #jalla fix. Burde velge en root finder som returnerer en skalar
    u_star = phi_l(p_star) #Finding velocity at intermediate state

    rho_l_star = rho_l*(1 + beta*p_star/p_l)/(p_star/p_l + beta) #Density at the left side of contact discontinuity
    rho_r_star = rho_r*(1 + beta*p_star/p_r)/(p_star/p_r + beta) #Density at the right side of contact discontinuity
    return u_star, rho_l_star, rho_r_star, p_star

    #Rarefaction structure
    #rarefactions in field 1:
def rarefaction1(xi,u_l,c_l,rho_l,p_l):
    gamma = 1.4
    u = ((gamma-1)*u_l + 2*(c_l+xi))/(gamma+1)
    rho = ((rho_l**gamma*(u - xi)**2)/(gamma*p_l))**(1/(gamma-1))
    p = (p_l/rho_l**gamma)*rho**gamma
    return rho,u,p

#rarefactions in field 3:
def rarefaction3(xi,u_r,c_r,rho_r,p_r):
    gamma = 1.4
    u = ((gamma-1)*u_r - 2*(c_r-xi))/(gamma+1)
    rho = ((rho_r**gamma*(u - xi)**2)/(gamma*p_r))**(1/(gamma-1))
    p = (p_r/rho_r**gamma)*rho**gamma
    return rho,u,p

def riemann_exact_solution(rho_l,u_l,p_l,rho_r,u_r,p_r, ni, L_x, t):

    gamma = 1.4
    beta = (gamma+1)/(gamma-1)
    c_l = math.sqrt(gamma*p_l/rho_l)
    c_r = math.sqrt(gamma*p_r/rho_r)

    u_star, rho_l_star, rho_r_star, p_star = riemann_exact(rho_l,u_l,p_l,rho_r,u_r,p_r)

    #Determine the wave speeds of the problem. In the case of a rarefaction, the speeds to the left and right of
    # the rarefaction is found
    if p_star <= p_l:
        #1-rarefaction. Note that this case is also able to handle the case when p_star == p_l, where the characteristic speed
        # becomes the same at both sides of the wave
        s1_l = u_l - c_l
        c_l_star = math.sqrt(gamma*p_star/rho_l_star)
        s1_r = u_star - c_l_star
    elif p_star > p_l:
        #1-shock
        s1 = (rho_l_star*u_star - rho_l*u_l)/(rho_l_star - rho_l)

    #contact discontinuity
    s2 = u_star

    if p_star <= p_r:
        #3-rarefaction. Note that this case is also able to handle the case when p_star == p_r, where the characteristic speed
        # becomes the same at both sides of the wave
        c_r_star = math.sqrt(gamma*p_star/rho_r_star)
        s3_l = u_star + c_r_star
        s3_r = u_r + c_r
    elif p_star > p_r:
        #3-shock
        s3 = (rho_r*u_r - rho_r_star*u_star)/(rho_r - rho_r_star)


    #Determine the solution as a function of x(t)
    dx = L_x/ni
    x = np.linspace(-(L_x/2 - dx/2),L_x/2-dx/2,ni) #centering the domain so that the shock happens at x = 0

    rho = np.ones(x.shape)*rho_r
    u = np.ones(x.shape)*u_r
    p = np.ones(x.shape)*p_r



    if p_star <= p_l: #1-rarefaction
        #left of 1-rarefaction
        ind = np.where(x <= s1_l*t)
        rho[ind] = rho_l
        u[ind] = u_l
        p[ind] = p_l
        #inside 1-rarefaction
        ind = np.where(np.logical_and(x > s1_l*t, x <= s1_r*t))
        xi = x[ind]/t
        rho_tmp,u_tmp,p_tmp = rarefaction1(xi,u_l,c_l,rho_l,p_l)
        rho[ind] = rho_tmp
        u[ind] = u_tmp
        p[ind] = p_tmp
        #left of contact discontinuity
        ind = np.where(np.logical_and(x > s1_r*t, x <= s2*t))
        rho[ind] = rho_l_star
        u[ind] = u_star
        p[ind] = p_star
    elif p_star > p_l: #1-shock
        #left of shock
        ind = np.where(x <= s1*t)
        rho[ind] = rho_l
        u[ind] = u_l
        p[ind] = p_l
        #right of shock
        ind = np.where(np.logical_and(x > s1*t, x <= s2*t))
        rho[ind] = rho_l_star
        u[ind] = u_star
        p[ind] = p_star

    if p_star <= p_r: #3-rarefaction
        #left of 3-rarefaction
        ind = np.where(np.logical_and(x > s2*t, x <= s3_l*t))
        rho[ind] = rho_r_star
        u[ind] = u_star
        p[ind] = p_star
        #inside 3-rarefaction
        ind = np.where(np.logical_and(x > s3_l*t, x <= s3_r*t))
        xi = x[ind]/t
        rho_tmp,u_tmp,p_tmp = rarefaction3(xi,u_r,c_r,rho_r,p_r)
        rho[ind] = rho_tmp
        u[ind] = u_tmp
        p[ind] = p_tmp
    elif p_star > p_r: #3-shock
        #left of 3-shock
        ind = np.where(np.logical_and(x > s2*t, x <= s3*t))
        rho[ind] = rho_r_star
        u[ind] = u_star
        p[ind] = p_star
        #right of shock is allready set
    return rho,u,p



def static_wall_exact(rho_l,u_l,p_l):

    #Calculates the state left of the wall, and assumes a left running shock wave
    gamma = 1.4
    beta = (gamma+1)/(gamma-1)
    c_l = math.sqrt(gamma*p_l/rho_l)
    phi = lambda p: u_l + 2*c_l/math.sqrt(2*gamma*(gamma-1))*(1 - p/p_l)/(math.sqrt(1 + beta*p/p_l))
    p_star = op.fsolve(phi,1)

    p_star = p_star[0]
    u_star = 0
    rho_l_star = (1+beta*p_star/p_l)/(p_star/p_l + beta)*rho_l
    s = (rho_l_star*0 - rho_l*u_l)/(rho_l_star - rho_l) #rankine hugoniot
    return rho_l_star, p_star, s


def piston_exact(rho_l,p_l,rho_r,p_r,u_b,ni,L_x,width_b,t):
    #this problem is solved by solving two separate riemann problems, where the velocity in the star region equals the
    #piston velocity u_b. piston is assumed to always move to the right, and the fluid is initially still, such that
    #the wave structure (1- rarefaction on left problem and 3-shock on right side can be assumed)
    dx = L_x/ni
    x = np.linspace(dx/2,L_x-dx/2,ni)
    g = 1.4
    beta = (g+1)/(g-1)
    u_l = 0
    u_r = 0
    u_star = u_b
    rho = np.ones(x.shape)*rho_l
    u = np.ones(x.shape)*u_l
    p = np.ones(x.shape)*p_l

    #left problem: 1-rarefaction

    c_l = math.sqrt(g*p_l/rho_l)
    phi_l = lambda p: u_l + 2*c_l/(g-1)*(1 - (p/p_l)**((g-1)/(2*g))) - u_star
    p_star = op.fsolve(phi_l,1)
    p_star = p_star[0]
    rho_star = (p_star/p_l)**(1/g)*rho_l
    c_star = math.sqrt(g*p_star/rho_star)
    s1_l = u_l - c_l
    s1_r = u_star - c_star

    x0 = 0.5*(L_x - width_b)
    ind = np.where(np.logical_and(x > x0 + s1_l*t, x <= x0 + s1_r*t))
    xi = (x[ind]-x0)/t
    rho_tmp,u_tmp,p_tmp = rarefaction1(xi,u_l,c_l,rho_l,p_l)
    rho[ind] = rho_tmp
    u[ind] = u_tmp
    p[ind] = p_tmp

    ind = np.where(np.logical_and(x > x0 + s1_r*t, x <= x0 + u_b*t))
    rho[ind] = rho_star
    u[ind] = u_star
    p[ind] = p_star

    #right problem: 3- shock
    x0 = 0.5*(L_x + width_b)
    c_r = math.sqrt(g*p_r/rho_r)
    phi_r = lambda p: u_r - 2*c_r/math.sqrt(2*g*(g-1))*(1 - p/p_r)/(math.sqrt(1 + beta*p/p_r)) - u_star
    p_star = op.fsolve(phi_r,p_r+1)
    p_star = p_star[0]
    rho_star_r = (1+beta*p_star/p_r)/(p_star/p_r+beta)*rho_r
    s3 = (rho_r*0 - rho_star_r*u_star)/(rho_r - rho_star_r)

    ind = np.where(np.logical_and(x >= x0+u_star*t, x<x0 + s3*t))
    rho[ind] = rho_star_r
    u[ind] = u_star
    p[ind] = p_star

    ind = np.where(x >= x0 + s3*t)
    rho[ind] = rho_r
    u[ind] = u_r
    p[ind] = p_r

    #inside piston
    x0 = 0.5*(L_x - width_b)
    ind = np.where(np.logical_and(x > x0+u_b*t, x <= x0+width_b+u_b*t))
    rho[ind] = float('nan')
    u[ind] = float('nan')
    p[ind] = float('nan')

    return rho,u,p




