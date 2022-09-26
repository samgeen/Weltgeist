"""
Manage gravity in the simulation volume
Sam Geen, March 2018
"""

import numpy as np

from . import integrator, units

# Use gravity?
gravity_on = False

# Add gravity to the whole volume?
# - if true, apply gravity to the whole volume
# - if false, apply gravity only to perturbed gas (to simulate turbulent support)
wholevol = True
centralmass = 0.0

def calculate_gravity():
    """
    Inject gravitational *acceleration* onto the grid
    NOTE: This is not fully tested yet, so be careful
    """
    hydro = integrator.Integrator().hydro
    # g = G*M/r^2
    ncells = hydro.ncells
    mass = hydro.mass[0:ncells]
    # Mass should be the cumulative mass of the grid *except* this cell
    mass = np.cumsum(mass) - mass
    mass += centralmass
    radius = hydro.x[0:ncells]
    # -ve because it points towards the origin from +ve radius
    grav = mass*0.0 # initialising to zero
    mask = radius > 0.0
    grav[mask] = -units.G * mass[mask] / (radius[mask])**2
    #grav[hydro.nH[0:ncells]<1e-6] = 0.0 # Prevent small timesteps
    if not wholevol:
        # Do this
        print("Implement wholevol=False for just the shocked region in gravity.py")
        raise NotImplementedError
    else:
        # Add gravitational force to the whole grid
        hydro.grav[0:ncells] = grav