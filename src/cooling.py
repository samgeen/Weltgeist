"""
Controls gas cooling
Sam Geen, March 2018
"""

import cooling_module, integrator

# Use cooling?
cooling_on = True

def solve_cooling(dt):
    """
    Solve the cooling step

    Parameters
    ----------

    dt : float
        timestep in seconds
    """
    # Initialise
    hydro = integrator.Integrator().hydro
    ncell = hydro.ncells
    nH = hydro.nH[0:ncell]
    T2 = hydro.T[0:ncell]
    xhii = hydro.xhii
    zsolar = hydro.Zsolar[0:ncell]
    gamma = hydro.gamma
    # Solve the change in temperature and modify the grid
    # Uses the model by Audit & Hennebelle (2005)
    # Used in the FRIG project by Patrick Hennebelle
    dT2 = cooling_module.solve_cooling_frig(nH,T2,zsolar,dt,gamma,ncell)
    T2 += dT2
    #hydro.T[xhii < 0.1] = T2[xhii < 0.1]
    hydro.T[0:ncell] = T2 * (1.0-xhii)
    # Make sure the resulting temperatures aren't negative
    CheckTemperature()

def CheckTemperature():
    """
    Check for and fix super low tempertaures
    """
    hydro = integrator.Integrator().hydro
    ncell = hydro.ncells
    T2 = hydro.T[0:ncell]
    # Extra check for low temperatures
    T2[T2 < 1.0] = 1.0
    hydro.T[0:ncell] = T2
