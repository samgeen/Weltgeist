"""
Spread radiation from a central source through the grid
Sam Geen, March 2018
"""

import numpy as np
from . import sources, integrator, units, ionisedtemperatures
from . import raytracing

# Dust cross section to use (Draine suggests 1e-21 cm^2 / H)
sigmaDust = 1e-21 

# Force photoionised temperature to Tion? WARNING: do not turn on if you have winds, useful for simple tests
forceTion = False

# Turn on radiation tracing?
radiation_on = True

def alpha_B_HII(temperature):
    """
    Calculate the HII recombination rate
    This is the rate at which ionised hydrogen recombines 
      into neutral hydrogen
    Total recombinations per second per unit volume = alpha * ne * nH
    ne = electron number density
    nH = hydrogen number density

    Parameters
    ----------

    temperature: float
        Temperature in K

    Returns
    -------

    alpha_B_HII : float
        The recombination rate
    """
    # HII recombination rate
    # input  : T in K
    # output : HII recombination rate (in cm3 / s)
    l = 315614./temperature
    a = 2.753e-14 * l**1.5 / (1. + (l/2.74)**0.407)**2.242
    return a             

def IonisedGasTemperature(Teff, metal):
    """
    Calculate the temperature of the ionised gas in 

    Teff : float
        Effective tempertaure of the star in K

    metalin : float
        Metallicity of the gas
    """
    Tion = ionisedtemperatures.FindTemperature(Teff, metal)
    return Tion

def trace_radiation(Lionising, Lnonionising, Eionising, Tion, doRadiationPressure):
    """
    Trace a ray through the spherical grid and ionise everything in the way
    Use very simple instant ionisation/recombination model
    d(nphotons)/dt = integrate(4*pi*r^2*nH(r)^2*alpha_B*dr)

    Parameters
    ----------
    Lionising, Lnonionising : float
        Luminosity of (non-)ionising photons emitted (erg/s)

    Eionising : float
        Average energy of the ionising photons (erg)

    Tion : float
        Equilibrium temperature of photoionised gas in K

    doRadiationPressure : bool
        Do we use radiation pressure?
    """

    global sigmaDust, forceTion

    # Note: we need to trace radiation even if there are no photons to calculate recombination

    hydro = integrator.Integrator().hydro
    dt = integrator.Integrator().dt
    # Photon emission rate
    QH = 0.0
    if Eionising > 0.0:
        QH = Lionising / Eionising

    # Find total recombinations from the centre outwards
    # Cool everything below gracefact*Tion to Tion to limit wiggles
    #gracefact = 1.001 
    gracefact = 2.0 # Use a larger gracefact to damp wiggles in temperature
    if Tion > 0.0:
        alpha_B = alpha_B_HII(Tion) 
    else:
        # Fallback to just using the gas temperature for recombining gas
        alpha_B = alpha_B_HII(hydro.T[:]) 
    alpha_B = 2.7e-13
    nx = hydro.ncells
    T = hydro.T[0:nx]
    x = hydro.x[0:nx]
    dx = hydro.dx
    c = units.c

    # Set up radiation tracing
    hydro.Qion[0] = QH
    if sigmaDust is not None:
        hydro.sigmaDust[0:nx] = sigmaDust # cm^2 / H based on Draine+ 2011
    hydro.sigmaDust[T > 1e5] = 0.0 # simplicity hack - remove dust from hot gas

    # Rate of recombinations per radial element
    drecombinationsdr = 4.0*np.pi*(x[0:nx]+dx)**2.0 * hydro.nH[0:nx]**2.0 * alpha_B
    drecombinationsdr[T > 1e5] = 0.0 # Assume collisionally ionised

    # Total number of ionising photon absorptions
    recombinations = np.cumsum(drecombinationsdr) * dx

    # Calculate new Qion along ray (function will modify hydro.Qion)
    raytracing.trace_radiation(dx,hydro.Qion[0:nx],hydro.sigmaDust[0:nx],hydro.nH[0:nx],drecombinationsdr,nx)

    # Calculate optical depth for non-ionising radiation
    opticalDepth = np.cumsum(hydro.nH[0:nx] * hydro.sigmaDust[0:nx]) * dx

    # Calculate recombinations first
    numatomspercell = hydro.nH*hydro.vol
    numionspercell = numatomspercell * hydro.xhii
    numionspercell -= recombinations*dt
    numionspercell[np.where(numionspercell < 0)] = 0.0
    hydro.xhii[0:nx] = numionspercell / numatomspercell

    # Calculate which gas should be ionised
    ionised = np.where(recombinations < QH)[0]

    # Note: Set xhii BEFORE T as this affects the final gas pressure (which takes into account free electron density)

    # Ionise all fully-ionised cells
    edge = 0
    if len(ionised) > 0:
        if not forceTion:
            toionise = (recombinations < QH)*(hydro.T[0:nx]<gracefact*Tion)
        else:
            # Force Tion to be set even if temperature is elevated
            toionise = (recombinations < QH)
        #print(len(hydro.xhii[ionised]) / len(hydro.xhii[toionise]))
        hydro.xhii[ionised] = 1.0
        hydro.T[toionise] = Tion
        edge = ionised[-1]+1

    # Ionise the partially ionised frontier cell provided it's not outside the box
    # Also check for case where there are no photons and the edge should not exist
    if edge < nx and QH > 0.0:
        Qextra = recombinations[edge] - QH
        if edge > 0:
            fracion = Qextra / (recombinations[edge]-recombinations[edge-1])
        else:
            fracion = Qextra / recombinations[edge]
        hydro.xhii[edge] = fracion
        # Fractionally heat the edge cell as if the ionisation front is sharp
        if hydro.T[edge] < Tion:
            hydro.T[edge] = hydro.T[edge]*(1.0-fracion) + fracion*Tion
    
    # Now do radiation pressure
    if doRadiationPressure:
        # See Draine (2011) equation 1 to give dP / dR
        mp = units.mH / units.X
        vol = hydro.vol[0:nx]

        # Calculate F = PA rather than P since there is a singularity in 1/(4 pi r^2) at r=0

        # Add dust contribution
        dFdrDust = hydro.nH[0:nx] * hydro.sigmaDust * \
            (Lnonionising * np.exp(-opticalDepth) + hydro.Qion[0:nx] * Eionising) / c

        # Add the direct radiation pressure contribution
        dFdrDirect = alpha_B * (hydro.nH[0:nx] * hydro.xhii[0:nx])**2 * Eionising / c * 4 * np.pi * x**2 

        dFdr = dFdrDust + dFdrDirect

        # Convert to outward acceleration from gravity
        # Note: we're using grav = acceleration here
        # P = F/A = ma / A
        # So a = P A / m
        # A = 4 pi r^2
        # m = rho * vol = nH * mp * vol
        # So correct Draine equation 1 by
        # 4 pi r^2 / (nH * vol * mp)

        # Convert pressure to acceleration
        hydro.grav[0:nx] = dFdr / hydro.mass[0:nx]

        # Thermal pressure handled by hydro solver
    # So we're done!
