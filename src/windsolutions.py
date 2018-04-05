'''
Composite Spitzer solution
Sam Geen, September 2016
'''

import numpy as np

import gravity, init, radiation, units

Myr = 3.1557e13 # seconds
pc = 3.0857e18 # cm
mH  = 1.66e-24 # g
mp = mH / units.X # g (mH / X)

def AdiabaticWind(lum,rho,time,model="Castor"):
    '''
    Adiabatic solution
    lum = luminosity in ergs/s
    rho = density in H/cc
    time = time in Myr
    model = which model to use ("Castor","Avedisova")
    '''
    # Convert units
    rho = rho * mp
    #time = time * Myr
    lum = lum # should be ok
    wconst = 0.88
    if model == "WeaverIntermediate":
        wconst = 0.76
    if model == "Avedisova":
        wconst = 1.02
    #wconst = 0.968
    # Equation 6, Castor et al 1975 (see also Avedisova 1972, Weaver 1977)
    return wconst*(lum * time ** 3 / rho)**0.2

def SpitzerSolution(Sphotons,n0,time):
    '''
    Spitzer solution for a photoionisation front
    '''
    Tion = 8400.0 # K, value used in Geen+ 2015b
    ci = np.sqrt(Tion * init.gamma * units.kB / mp)
    alpha_B = radiation.alpha_B_HII(Tion)
    rs = (Sphotons / (4.0/3.0 * np.pi * alpha_B * n0**2))**(1.0/3.0)
    # Hosokawa & Inutsuka 2004 approx
    #rspitzer = rs  * (1.0 + 7.0/4.0 * np.sqrt(4.0/3.0) * ci / rs * time)**(4.0/7.0)
    rspitzer = rs  * (1.0 + 7.0/4.0 * ci / rs * time)**(4.0/7.0)
    return rspitzer

def CollapseSolution(rho,rho0):
    '''
    Calculate a free-fall collapse solution
    rho - density to calculate time at
    rho0 - initial density (used to calculate tff)
    Sam Geen, March 2018
    '''
    tff = np.sqrt(3.0*np.pi / (32.0 * units.G * rho0))
    X = (rho / rho0)**(-0.5)
    #print X, tff, rho, rho0
    t = tff * 2.0 / np.pi * (np.arccos(np.sqrt(X)) + np.sqrt(X * (1.0-X)))
    return t

def CollapseSolutionPosition(x,x0):
    '''
    Calculate a free-fall collapse solution
    x - position to calculate time at
    x0 - initial position
    Sam Geen, March 2018
    '''
    X = x/x0
    t = (np.arccos(np.sqrt(X)) + np.sqrt(X * (1.0-X))) * x0**1.5 / np.sqrt(2.0*units.G*gravity.centralmass)
    return t