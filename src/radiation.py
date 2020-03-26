"""
Spread radiation from a central source through the grid
Sam Geen, March 2018
"""

import numpy as np

import sources, integrator, units, windsolutions

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

def trace_radiation(totalphotons, Tion=8400.0):
    """
    Trace a ray through the spherical grid and ionise everything in the way
    Use very simple instant ionisation/recombination model
    d(nphotons)/dt = integrate(4*pi*r^2*nH(r)^2*alpha_B*dr)

    Parameters
    ----------
    totalphotons : float
        Number of ionising photons emitted

    Tion : float
        Equilibrium temperature of photoionised gas in K
        default: 8400 K, value used in Geen+ 2015b
    """

    hydro = integrator.Integrator().hydro
    dt = integrator.Integrator().dt
    # Photon emission rate
    QH = totalphotons / dt

    # No photons? Don't bother
    if QH == 0:
        return

    #Dust Parameters
    sigma_dust = 1.0
    dgr = 1.0
    # Find total recombinations from the centre outwards
    Sphotons_dust = 0
    gracefact = 2.0 # Cool everything below 2*Tion to Tion to limit wiggles
    alpha_B = alpha_B_HII(Tion) 
    nx = hydro.ncells
    
    recombinations = np.cumsum(4.0*np.pi*(hydro.x[0:nx]+hydro.dx)**2.0 * hydro.nH[0:nx]**2.0 * alpha_B * hydro.dx)

    dTau_dust = hydro.nH[0:nx] * sigma_dust * dgr * hydro.dx
    # This is wrong. All of Sphotons is not available to dust. 
    # Need to calculate effective cross section for each species and
    # do both at the same time. 
    dN_dust = Sphotons_dust*(1-np.exp(-dTau_dust))

    ionised = np.where(recombinations + dN_dust*0 < QH)[0]
    # Ionise all fully-ionised cells
    edge = 0
    if len(ionised) > 0:
        toionise = (recombinations < QH)*(hydro.T[0:nx]<gracefact*Tion)
        hydro.xhii[ionised] = 1.0
        hydro.T[toionise] = Tion
        edge = ionised[-1]+1
    # Ionise the partially ionised frontier cell
    Sextra = recombinations[edge] - QH
    if edge > 0:
        fracion = Sextra / (recombinations[edge]-recombinations[edge-1])
    else:
        fracion = Sextra / recombinations[edge]
    hydro.xhii[edge] = fracion
    # Fractionally heat the edge cell as if the ionisation front is sharp
    if hydro.T[edge] < Tion:
        hydro.T[edge] = hydro.T[edge]*(1.0-fracion) + fracion*Tion
    # Done!

if __name__=="__main__":
    init.init()
    trace_radiation(1e48)