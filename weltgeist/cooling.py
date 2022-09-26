"""
Controls gas cooling
Sam Geen, March 2018
"""

from . import cooling_module, integrator

import numpy as np

# Use cooling?
cooling_on = False

# Mask out the contact discontinuity to cut out numerical cooling in this region
# Depends on your opinion of how this subgrid physics works
maskContactDiscontinuity = False

def TemperatureChange(dt):
    """
    Calculate the temperature change needed for each cell

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
    # Solve the change in temperature
    # Uses the model by Audit & Hennebelle (2005)
    # Used in the FRIGG project by Patrick Hennebelle
    dT2 = cooling_module.solve_cooling_frig(nH,T2,zsolar,dt,gamma,ncell)
    # Remove cooling blip in very centre that might be eating energy
    dT2[0:1] = 0.0
    return dT2

def MaskContactDiscontinuityV1(dT2):
    """
    Find the location of contact discontinuity to remove cooling from
    Prevents artificial mixing by numerical diffusion
    Version 1: just grabs the cells around the CD
    
    Parameters
    ----------

    dT2 : array
        temperature change to mask out
    """
    hydro = integrator.Integrator().hydro
    ncell = hydro.ncells
    shock = np.where(hydro.T[0:ncell] >= 1e6)[0]
    if len(shock) > 0:
        edge = shock[-1]
        dT2[edge-2:edge+3] = 0.0
    return dT2

def MaskContactDiscontinuityV2(dT2):
    """
    Find the location of contact discontinuity to remove cooling from
    Prevents artificial mixing by numerical diffusion

    Version 2: identifies whole numerically-smeared CD and tries to reconstruct it
    
    Parameters
    ----------

    dT2 : array
        temperature change to mask out
    Tion : float
        ionised gas temperature
    """
    hydro = integrator.Integrator().hydro
    ncell = hydro.ncells
    T = hydro.T[0:ncell]
    shock = np.where(hydro.T[0:ncell] >= 1e6)[0]
    if len(shock) > 0:
        inneredge = shock[-1]
        # Look for end of CD
        iCD = inneredge
        while ((T[iCD] - T[iCD+1]) > 1.0 and T[iCD] > 1e4):
            iCD += 1
        outeredge = iCD
        # Look for start of CD
        iCD = inneredge
        while ((T[iCD-1] - T[iCD] and inneredge >= 1) > 1.0):
            iCD -= 1
        inneredge = max(iCD-1,0)
        #print(inneredge, outeredge, dT2[inneredge:outeredge+2])
        
        #print(inneredge, outeredge, hydro.PMagnetic[inneredge:outeredge+2] / hydro.PThermal[inneredge:outeredge+2])
        dT2[inneredge:outeredge+1] = dT2[inneredge]
    return dT2


def solve_cooling(dt):
    """
    Solve the cooling step

    Parameters
    ----------

    dt : float
        timestep in seconds
    """
    # Solve the change in temperature and modify the grid
    hydro = integrator.Integrator().hydro
    ncell = hydro.ncells
    dT2 = TemperatureChange(dt)
    # Mask wind shock to prevent numerical diffusion cooling effects
    if maskContactDiscontinuity:
        dT2 = MaskContactDiscontinuityV2(dT2)
    hydro.T[0:ncell] += dT2
    # Note that the radiation module, if it runs, will heat the 
    #  photoionised gas back up
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
