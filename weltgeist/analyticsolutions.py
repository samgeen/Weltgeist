"""
Analytic solutions for test problems
Sam Geen, September 2016
"""

import numpy as np
import scipy.interpolate

from . import gravity, radiation, units, integrator

Myr = 3.1557e13 # seconds
yr = 1e-6 * Myr
pc = 3.0857e18 # cm
mH  = 1.66e-24 # g

def SedovTaylorSolution(time,energy,density):
    """
    Adiabatic solution for a blastwave
    See http://www.astronomy.ohio-state.edu/~ryden/ast825/ch5-6.pdf

    Parameters
    ----------

    time : float/array
        time(s) to calculate at in seconds
    energy : float
        energy of the blastwave in erg
    density : float
        density of the background in g cm^-3

    Returns
    -------

    radius : float/array
        a radius for each time inputted
    """ 
    return 1.17 * (energy * time**2 / density)**0.2


def AdiabaticWind(lum,ml,rho,time,model="Castor"):
    """
    Adiabatic solution for a constant wind
    Includes various models for how the wind energy behaves
    
    Parameters
    ----------

    lum = luminosity in ergs/s
    ml = mass loss rate in g/s
    rho = density in H/cc
    time = time in s
    model = which model to use ("Castor","Avedisova")
    """
    # Convert units
    rho = rho * units.mp
    #time = time * Myr
    lum = lum # should be ok
    mom = np.sqrt(2.0*lum*ml)
    # Castor/Weaver, no cooling
    wconst = 0.88
    tcool = yr * 3e2 # HACK FIT TO IFFRIG !!!! FIX THIS!!!!
    if model == "WeaverIntermediate":
        # Weaver cooled model
        wconst = 0.76
    elif model == "Avedisova":
        # Avedisova, cooling
        wconst = 1.02
    elif model == "CapriottiAdiabatic":
        # Capriotti/Weaver, no cooling
        wconst = 0.86
    elif model == "CooledNoAcc":
        # Different solution, momentum-dominated, assume no acceleration
        fcool = np.sqrt(2.0*mom/(np.pi*rho))
        rcool = 0.76 * (lum * tcool ** 3 / rho)**0.2
        return np.sqrt(rcool*rcool + fcool * (time - tcool))
    elif model == "Cooled":
        # Different solution, momentum-dominated, with acceleration
        fcool = np.sqrt(2.0*mom/(np.pi*rho))
        # Set rcool, vcool for Weaver cooled model
        rcool = 0.76 * (lum * tcool ** 3 / rho)**0.2
        vcool = 0.6 * rcool / tcool # dr/dt = 3/5 r/t since r \propto t^{3/5}
        tnew = time - tcool
        r4 = rcool**4 + (3.0 * mom) / (2.0 * np.pi * rho) * tnew**2.0 + 4.0 * rcool**3.0 * vcool * tnew
        if r4 < 0.0:
            r4 = 0.0
        #print rcool
        #import pdb; pdb.set_trace()
        return r4**0.25
    elif model == "SlowCool":
        # Assume bubble retains L2 * tcool energy, Sedov-like
        rcool = 0.76 * (lum * tcool ** 3 / rho)**0.2
        if time > tcool:
            return 0.7 *(lum * tcool * (time)**2 / rho)**0.2# + rcool
        else:
            return 0.0
    #wconst = 0.968
    # Equation 6, Castor et al 1975 (see also Avedisova 1972, Weaver 1977)
    return wconst*(lum * time ** 3 / rho)**0.2

def SpitzerSolution(QH,n0,time, Tion = 1e4):
    """
    Spitzer solution for a photoionisation front
    From the Spitzer 1978 book on the interstellar medium

    Parameters
    ----------
    QH : float
        photon emission rate in photons/second
    n0 : float
        hydrogen number density in cm^-3
    time : float or numpy array
        time in s
    Tion : float
        ionised gas temperature in K

    Returns
    -------
    radius : float or numpy array
        radius in cm (same length as time)
    """
    ci = np.sqrt(2.0 * Tion * units.kB / units.mp)
    alpha_B = radiation.alpha_B_HII(Tion)
    rs = (QH / (4.0/3.0 * np.pi * alpha_B * n0**2))**(1.0/3.0)
    rspitzer = rs  * (1.0 + 7.0/4.0 * ci / rs * time)**(4.0/7.0)
    return rspitzer

def SpitzerDensity(QH,n0,time,Tion = 8400.0):
    """
    Spitzer solution for the ionised gas density
    Derived from the solution above

    Parameters
    ----------
    QH : float
        photon emission rate in photons/second
    n0 : float
        hydrogen number density in cm^-3
    time : float or numpy array
        time in s

    Returns
    -------
    ni : float or numpy array
        hydrogen number density in ionised gas
    """
    # Get the radius
    ri = SpitzerSolution(QH,n0,time,Tion)
    # Get density from photoionisation equilibrium
    alpha_B = radiation.alpha_B_HII(Tion)
    ni = np.sqrt(3.0 * QH / (4.0 * np.pi * alpha_B * ri**3.0))
    return ni

def HosokawaInutsuka(QH,n0,time, Tion = 1e4):
    """
    Hosokawa & Inutsuka solution for a photoionisation front
    From Hosokawa & Inutsuka (2006)
    Used by the Starbench paper (Bisbas et al. 2015)

    Parameters
    ----------
    QH : float
        photon emission rate in photons/second
    n0 : float
        hydrogen number density in cm^-3
    time : float or numpy array
        time in s
    Tion : float
        ionised gas temperature in K

    Returns
    -------
    radius : float or numpy array
        radius in cm (same length as time)
    """
    ci = np.sqrt(2.0 * Tion * units.kB / units.mp)
    alpha_B = radiation.alpha_B_HII(Tion)
    rs = (QH / (4.0/3.0 * np.pi * alpha_B * n0**2))**(1.0/3.0)
    # Hosokawa & Inutsuka 2006 approx
    rhi = rs  * (1.0 + 7.0/4.0 * np.sqrt(4.0/3.0) * ci / rs * time)**(4.0/7.0)
    return rhi

def StarbenchExtracted(QH,n0,time, Tion = 1e4):
    """
    The "Starbench" semi-analytic solution for late phase evolution in their test
    From Bisbas et al. (2015), equation 28, Figure 5

    Parameters
    ----------
    QH : float
        photon emission rate in photons/second
    n0 : float
        hydrogen number density in cm^-3
    time : float or numpy array
        time in s
    Tion : float
        ionised gas temperature in K

    Returns
    -------
    radius : float or numpy array
        radius in cm (same length as time)
    """
    # Solving this was too hard so instead I'm just reading scraped data from the figure
    '''
    ci = np.sqrt(2.0 * Tion * units.kB / units.mp)
    alpha_B = radiation.alpha_B_HII(Tion)
    rs = (QH / (4.0/3.0 * np.pi * alpha_B * n0**2))**(1.0/3.0)
    RI = 0
    RII = 0
    fSB = 1.0 - 0.733 * np.exp(-time / wunits.Myr)
    # Starbench equation
    rSB = RII + fSB * (RI - RII)
    '''
    # Even values are x coordinates
    traw = np.array(starbench_late_raw[0::2]) * units.Myr
    rraw = np.array(starbench_late_raw[1::2]) * units.pc
    # Make an interpolation function to remap to the requested times
    rfunc = scipy.interpolate.interp1d(traw,rraw,fill_value="extrapolate")
    # Do that thing I said
    rSB = rfunc(time)
    # And return it
    return rSB

def CollapseSolution(rho,rho0):
    """
    Calculate a free-fall collapse solution
    rho - density to calculate time at in g/cm^3
    rho0 - initial density (used to calculate tff) in g/cm^3
    Sam Geen, March 2018
    """
    tff = np.sqrt(3.0*np.pi / (32.0 * units.G * rho0))
    X = (rho / rho0)**(-0.5)
    #print X, tff, rho, rho0
    t = tff * 2.0 / np.pi * (np.arccos(np.sqrt(X)) + np.sqrt(X * (1.0-X)))
    return t

def CollapseSolutionPosition(x,x0):
    """
    Calculate a free-fall collapse solution
    x - position to calculate time at in cm
    x0 - initial position in cm
    Taken from the ultra-scientific source:
    https://en.wikipedia.org/wiki/Equations_for_a_falling_body
    Sam Geen, March 2018
    """
    X = x/x0
    t = (np.arccos(np.sqrt(X)) + np.sqrt(X * (1.0-X))) * x0**1.5 / np.sqrt(2.0*units.G*gravity.centralmass)
    return t


'''
Starbench "Late" solution data dump
Extracted by drawing over Bisbas+ (2015) Figure 5
I could have solved the equations the long way or asked Thomas, but oh well
'''

starbench_late_raw = \
[0.0, 0.3195719952461076,
0.044776119402984926, 0.7391030884872443,
0.12686567164179097, 1.1873196067233636,
0.1865671641791044, 1.4185316951951705,
0.2499999999999999, 1.606283474046518,
0.3171641791044775, 1.77953047585238,
0.38432835820895517, 1.923821945083267,
0.4626865671641791, 2.0535546157526703,
0.5373134328358209, 2.1688365308925857,
0.6119402985074627, 2.2406851471700397,
0.6902985074626866, 2.3125067526894956,
0.7686567164179104, 2.3698505919214634,
0.8470149253731344, 2.4271944311534313,
0.9328358208955226, 2.4555287162944284,
1.0074626865671643, 2.4839440337094203,
1.0932835820895526, 2.4978005525629294,
1.1716417910447765, 2.511711092932437,
1.25, 2.5256216333019443,
1.3358208955223883, 2.5250003858679664,
1.4141791044776122, 2.524433159949986,
1.4962686567164178, 2.523838923274007,
1.5746268656716418, 2.5232716973560274,
1.6567164179104479, 2.5081996943925615,
1.7350746268656718, 2.507632468474581,
1.8134328358208958, 2.4925874762691143,
1.8955223880597019, 2.4919932395931355,
1.977611940298508, 2.4841601197734136,
2.057835820895523, 2.47634050533269,
2.139925373134329, 2.468507385512968,
2.2164179104477615, 2.4607147818302435,
2.300373134328358, 2.4456292734877785,
2.378731343283582, 2.437823164426055,
2.457089552238806, 2.430017055364331,
2.5410447761194033, 2.414931547021866,
2.6212686567164187, 2.4071119325811416,
2.701492537313433, 2.3992923181404198,
2.7817164179104488, 2.391472703699697,
2.8619402985074642, 2.3908919724027164,
2.942164179104478, 2.3758334748182515,
2.992537313432836, 2.375468829585264]
