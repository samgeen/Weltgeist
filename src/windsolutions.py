'''
Composite Spitzer solution
Sam Geen, September 2016
'''

import numpy as np

alpha_B = 2e-13 # idk
Myr = 3.1557e13 # seconds
pc = 3.0857e18 # cm
mH  = 1.66e-24 # g
mp = mH / 0.76 # g (mH / X)

S0 = 1e48
S1 = 1e51
t1 = 1.0

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
    wconst = 0.76
    if model == "Avedisova":
        wconst = 1.02
    # Equation 6, Castor et al 1975 (see also Avedisova 1972, Weaver 1977)
    return wconst*(lum * time ** 3 / rho)**0.2


